#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <ctype.h>
#include "sam.h"
#include "faidx.h"
#include "bam_maqcns.h"
#include "khash.h"
#include "glf.h"
#include "kstring.h"

typedef int *indel_list_t;
KHASH_MAP_INIT_INT64(64, indel_list_t)

#define BAM_PLF_SIMPLE     0x01
#define BAM_PLF_CNS        0x02
#define BAM_PLF_INDEL_ONLY 0x04
#define BAM_PLF_GLF        0x08
#define BAM_PLF_VAR_ONLY   0x10
#define BAM_PLF_2ND        0x20
#define BAM_PLF_RANBASE    0x40
#define BAM_PLF_1STBASE    0x80
#define BAM_PLF_ALLBASE    0x100
#define BAM_PLF_READPOS    0x200

typedef struct {
	bam_header_t *h;
	bam_maqcns_t *c;
	bam_maqindel_opt_t *ido;
	faidx_t *fai;
	khash_t(64) *hash;
	uint32_t format;
	int tid, len, last_pos;
	int mask;
    int max_depth;  // for indel calling, ignore reads with the depth too high. 0 for unlimited
	char *ref;
	glfFile fp_glf; // for glf output only
} pu_data_t;

char **__bam_get_lines(const char *fn, int *_n);
void bam_init_header_hash(bam_header_t *header);
int32_t bam_get_tid(const bam_header_t *header, const char *seq_name);

static khash_t(64) *load_pos(const char *fn, bam_header_t *h)
{
	char **list;
	int i, j, n, *fields, max_fields;
	khash_t(64) *hash;
	bam_init_header_hash(h);
	list = __bam_get_lines(fn, &n);
	hash = kh_init(64);
	max_fields = 0; fields = 0;
	for (i = 0; i < n; ++i) {
		char *str = list[i];
		int chr, n_fields, ret;
		khint_t k;
		uint64_t x;
		n_fields = ksplit_core(str, 0, &max_fields, &fields);
		if (n_fields < 2) continue;
		chr = bam_get_tid(h, str + fields[0]);
		if (chr < 0) {
			fprintf(stderr, "[load_pos] unknown reference sequence name: %s\n", str + fields[0]);
			continue;
		}
		x = (uint64_t)chr << 32 | (atoi(str + fields[1]) - 1);
		k = kh_put(64, hash, x, &ret);
		if (ret == 0) {
			fprintf(stderr, "[load_pos] position %s:%s has been loaded.\n", str+fields[0], str+fields[1]);
			continue;
		}
		kh_val(hash, k) = 0;
		if (n_fields > 2) {
			// count
			for (j = 2; j < n_fields; ++j) {
				char *s = str + fields[j];
				if ((*s != '+' && *s != '-') || !isdigit(s[1])) break;
 			}
			if (j > 2) { // update kh_val()
				int *q, y, z;
				q = kh_val(hash, k) = (int*)calloc(j - 1, sizeof(int));
				q[0] = j - 2; z = j; y = 1;
				for (j = 2; j < z; ++j)
					q[y++] = atoi(str + fields[j]);
			}
		}
		free(str);
	}
	free(list); free(fields);
	return hash;
}

// an analogy to pileup_func() below
static int glt3_func(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pu, void *data)
{
	pu_data_t *d = (pu_data_t*)data;
	bam_maqindel_ret_t *r = 0;
	int rb, *proposed_indels = 0;
	glf1_t *g;
	glf3_t *g3;

	if (d->fai == 0) {
		fprintf(stderr, "[glt3_func] reference sequence is required for generating GLT. Abort!\n");
		exit(1);
	}
	if (d->hash) { // only output a list of sites
		khint_t k = kh_get(64, d->hash, (uint64_t)tid<<32|pos);
		if (k == kh_end(d->hash)) return 0;
		proposed_indels = kh_val(d->hash, k);
	}
	g3 = glf3_init1();
	if (d->fai && (int)tid != d->tid) {
		if (d->ref) { // then write the end mark
			g3->rtype = GLF3_RTYPE_END;
			glf3_write1(d->fp_glf, g3);
		}
		glf3_ref_write(d->fp_glf, d->h->target_name[tid], d->h->target_len[tid]); // write reference
		free(d->ref);
		d->ref = fai_fetch(d->fai, d->h->target_name[tid], &d->len);
		d->tid = tid;
		d->last_pos = 0;
	}
	rb = (d->ref && (int)pos < d->len)? d->ref[pos] : 'N';
	g = bam_maqcns_glfgen(n, pu, bam_nt16_table[rb], d->c);
	memcpy(g3, g, sizeof(glf1_t));
	g3->rtype = GLF3_RTYPE_SUB;
	g3->offset = pos - d->last_pos;
	d->last_pos = pos;
	glf3_write1(d->fp_glf, g3);
    if (pos < d->len) {
        int m = (!d->max_depth || d->max_depth>n) ? n : d->max_depth;
		if (proposed_indels)
			r = bam_maqindel(m, pos, d->ido, pu, d->ref, proposed_indels[0], proposed_indels+1);
		else r = bam_maqindel(m, pos, d->ido, pu, d->ref, 0, 0);
	}
	if (r) { // then write indel line
		int het = 3 * n, min;
		min = het;
		if (min > r->gl[0]) min = r->gl[0];
		if (min > r->gl[1]) min = r->gl[1];
		g3->ref_base = 0;
		g3->rtype = GLF3_RTYPE_INDEL;
		memset(g3->lk, 0, 10);
		g3->lk[0] = r->gl[0] - min < 255? r->gl[0] - min : 255;
		g3->lk[1] = r->gl[1] - min < 255? r->gl[1] - min : 255;
		g3->lk[2] = het - min < 255? het - min : 255;
		g3->offset = 0;
		g3->indel_len[0] = r->indel1;
		g3->indel_len[1] = r->indel2;
		g3->min_lk = min < 255? min : 255;
		g3->max_len = (abs(r->indel1) > abs(r->indel2)? abs(r->indel1) : abs(r->indel2)) + 1;
		g3->indel_seq[0] = strdup(r->s[0]+1);
		g3->indel_seq[1] = strdup(r->s[1]+1);
		glf3_write1(d->fp_glf, g3);
		bam_maqindel_ret_destroy(r);
	}
	free(g);
	glf3_destroy1(g3);
	return 0;
}

static void pileup_seq(const bam_pileup1_t *p, int pos, int ref_len, const char *ref)
{
	if (p->is_head) printf("^%c", p->b->core.qual > 93? 126 : p->b->core.qual + 33);
	if (!p->is_del) {
		int j, rb, c = bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b), p->qpos)];
		rb = (ref && pos < ref_len)? ref[pos] : 'N';
		if (c == '=' || toupper(c) == toupper(rb)) c = bam1_strand(p->b)? ',' : '.';
		else c = bam1_strand(p->b)? tolower(c) : toupper(c);
		putchar(c);
		if (p->indel > 0) {
			printf("+%d", p->indel);
			for (j = 1; j <= p->indel; ++j) {
				c = bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b), p->qpos + j)];
				putchar(bam1_strand(p->b)? tolower(c) : toupper(c));
			}
		} else if (p->indel < 0) {
			printf("%d", p->indel);
			for (j = 1; j <= -p->indel; ++j) {
				c = (ref && (int)pos+j < ref_len)? ref[pos+j] : 'N';
				putchar(bam1_strand(p->b)? tolower(c) : toupper(c));
			}
		}
	} else putchar('*');
	if (p->is_tail) putchar('$');
}

static int pileup_func(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pu, void *data)
{
	pu_data_t *d = (pu_data_t*)data;
	bam_maqindel_ret_t *r = 0;
	int i, rb, rms_mapq = -1, *proposed_indels = 0;
	uint64_t rms_aux;
	uint32_t cns = 0;

	// if GLF is required, suppress -c completely
	if (d->format & BAM_PLF_GLF) return glt3_func(tid, pos, n, pu, data);
	// if d->hash is initialized, only output the sites in the hash table
	if (d->hash) {
		khint_t k = kh_get(64, d->hash, (uint64_t)tid<<32|pos);
		if (k == kh_end(d->hash)) return 0;
		proposed_indels = kh_val(d->hash, k);
	}
	// update d->ref if necessary
	if (d->fai && (int)tid != d->tid) {
		free(d->ref);
		d->ref = faidx_fetch_seq(d->fai, d->h->target_name[tid], 0, 0x7fffffff, &d->len);
		d->tid = tid;
	}
	rb = (d->ref && (int)pos < d->len)? d->ref[pos] : 'N';
	// when the indel-only mode is asked for, return if no reads mapped with indels
	if (d->format & BAM_PLF_INDEL_ONLY) {
		for (i = 0; i < n; ++i)
			if (pu[i].indel != 0) break;
		if (i == n) return 0;
	}
	// call the consensus and indel
	if (d->format & BAM_PLF_CNS) { // call consensus
		if (d->format & (BAM_PLF_RANBASE|BAM_PLF_1STBASE)) { // use a random base or the 1st base as the consensus call
			const bam_pileup1_t *p = (d->format & BAM_PLF_1STBASE)? pu : pu + (int)(drand48() * n);
			int q = bam1_qual(p->b)[p->qpos];
			int mapQ = p->b->core.qual < d->c->cap_mapQ? p->b->core.qual : d->c->cap_mapQ;
			uint32_t b = bam1_seqi(bam1_seq(p->b), p->qpos);
			cns = b<<28 | 0xf<<24 | mapQ<<16 | q<<8;
		} else if (d->format & BAM_PLF_ALLBASE) { // collapse all bases
			uint64_t rmsQ = 0;
			uint32_t b = 0;
			for (i = 0; i < n; ++i) {
				const bam_pileup1_t *p = pu + i;
				int q = p->b->core.qual < d->c->cap_mapQ? p->b->core.qual : d->c->cap_mapQ;
				b |= bam1_seqi(bam1_seq(p->b), p->qpos);
				rmsQ += q * q;
			}
			rmsQ = (uint64_t)(sqrt((double)rmsQ / n) + .499);
			cns = b<<28 | 0xf<<24 | rmsQ<<16 | 60<<8;
		} else cns = bam_maqcns_call(n, pu, d->c);
	}
    if ((d->format & (BAM_PLF_CNS|BAM_PLF_INDEL_ONLY)) && d->ref && pos < d->len) { // call indels
        int m = (!d->max_depth || d->max_depth>n) ? n : d->max_depth;
        if (proposed_indels) // the first element gives the size of the array
            r = bam_maqindel(m, pos, d->ido, pu, d->ref, proposed_indels[0], proposed_indels+1);
        else r = bam_maqindel(m, pos, d->ido, pu, d->ref, 0, 0);
	}
	// when only variant sites are asked for, test if the site is a variant
	if ((d->format & BAM_PLF_CNS) && (d->format & BAM_PLF_VAR_ONLY)) {
		if (!(bam_nt16_table[rb] != 15 && cns>>28 != bam_nt16_table[rb])) { // not a SNP
			if (!(r && (r->gt == 2 || strcmp(r->s[r->gt], "*")))) { // not an indel
				if (r) bam_maqindel_ret_destroy(r);
				return 0;
			}
		}
	}
	// print the first 3 columns
	printf("%s\t%d\t%c\t", d->h->target_name[tid], pos + 1, rb);
	// print consensus information if required
	if (d->format & BAM_PLF_CNS) {
		int ref_q, rb4 = bam_nt16_table[rb];
		ref_q = 0;
		if (rb4 != 15 && cns>>28 != 15 && cns>>28 != rb4) { // a SNP
			ref_q = ((cns>>24&0xf) == rb4)? cns>>8&0xff : (cns>>8&0xff) + (cns&0xff);
			if (ref_q > 255) ref_q = 255;
		}
		rms_mapq = cns>>16&0xff;
		printf("%c\t%d\t%d\t%d\t", bam_nt16_rev_table[cns>>28], cns>>8&0xff, ref_q, rms_mapq);
	}
	// print pileup sequences
	printf("%d\t", n);
	rms_aux = 0; // we need to recalculate rms_mapq when -c is not flagged on the command line
	for (i = 0; i < n; ++i) {
		const bam_pileup1_t *p = pu + i;
		int tmp = p->b->core.qual < d->c->cap_mapQ? p->b->core.qual : d->c->cap_mapQ;
		rms_aux += tmp * tmp;
		pileup_seq(p, pos, d->len, d->ref);
	}
	// finalize rms_mapq
	rms_aux = (uint64_t)(sqrt((double)rms_aux / n) + .499);
	if (rms_mapq < 0) rms_mapq = rms_aux;
	putchar('\t');
	// print quality
	for (i = 0; i < n; ++i) {
		const bam_pileup1_t *p = pu + i;
		int c = bam1_qual(p->b)[p->qpos] + 33;
		if (c > 126) c = 126;
		putchar(c);
	}
	if (d->format & BAM_PLF_2ND) { // print 2nd calls and qualities
		const unsigned char *q;
		putchar('\t');
		for (i = 0; i < n; ++i) {
			const bam_pileup1_t *p = pu + i;
			q = bam_aux_get(p->b, "E2");
			putchar(q? q[p->qpos + 1] : 'N');
		}
		putchar('\t');
		for (i = 0; i < n; ++i) {
			const bam_pileup1_t *p = pu + i;
			q = bam_aux_get(p->b, "U2");
			putchar(q? q[p->qpos + 1] : '!');
		}
	}
	// print mapping quality if -s is flagged on the command line
	if (d->format & BAM_PLF_SIMPLE) {
		putchar('\t');
		for (i = 0; i < n; ++i) {
			int c = pu[i].b->core.qual + 33;
			if (c > 126) c = 126;
			putchar(c);
		}
	}
	// print read position
	if (d->format & BAM_PLF_READPOS) {
		putchar('\t');
		for (i = 0; i < n; ++i) {
			int x = pu[i].qpos;
			int l = pu[i].b->core.l_qseq;
			printf("%d,", x < l/2? x+1 : -((l-1)-x+1));
		}
	}
	putchar('\n');
	// print the indel line if r has been calculated. This only happens if:
	// a) -c or -i are flagged, AND b) the reference sequence is available
	if (r) {
		printf("%s\t%d\t*\t", d->h->target_name[tid], pos + 1);
		if (r->gt < 2) printf("%s/%s\t", r->s[r->gt], r->s[r->gt]);
		else printf("%s/%s\t", r->s[0], r->s[1]);
		printf("%d\t%d\t", r->q_cns, r->q_ref);
		printf("%d\t%d\t", rms_mapq, n);
		printf("%s\t%s\t", r->s[0], r->s[1]);
		//printf("%d\t%d\t", r->gl[0], r->gl[1]);
		printf("%d\t%d\t%d\t", r->cnt1, r->cnt2, r->cnt_anti);
		printf("%d\t%d\n", r->cnt_ref, r->cnt_ambi);
		bam_maqindel_ret_destroy(r);
	}
	return 0;
}

int bam_pileup(int argc, char *argv[])
{
	int c, is_SAM = 0;
	char *fn_list = 0, *fn_fa = 0, *fn_pos = 0;
	pu_data_t *d = (pu_data_t*)calloc(1, sizeof(pu_data_t));
    d->max_depth = 0;
	d->tid = -1; d->mask = BAM_DEF_MASK;
	d->c = bam_maqcns_init();
	d->c->is_soap = 1; // change the default model
	d->ido = bam_maqindel_opt_init();
	while ((c = getopt(argc, argv, "st:f:cT:N:r:l:d:im:gI:G:vM:S2aR:PA")) >= 0) {
		switch (c) {
		case 'a': d->c->is_soap = 1; break;
		case 'A': d->c->is_soap = 0; break;
		case 's': d->format |= BAM_PLF_SIMPLE; break;
		case 't': fn_list = strdup(optarg); break;
		case 'l': fn_pos = strdup(optarg); break;
		case 'f': fn_fa = strdup(optarg); break;
		case 'T': d->c->theta = atof(optarg); break;
		case 'N': d->c->n_hap = atoi(optarg); break;
		case 'r': d->c->het_rate = atof(optarg); d->ido->r_snp = d->c->het_rate; break;
		case 'M': d->c->cap_mapQ = atoi(optarg); break;
		case 'd': d->max_depth = atoi(optarg); break;
		case 'c': d->format |= BAM_PLF_CNS; break;
		case 'i': d->format |= BAM_PLF_INDEL_ONLY; break;
		case 'v': d->format |= BAM_PLF_VAR_ONLY; break;
		case 'm': d->mask = strtol(optarg, 0, 0); break;
		case 'g': d->format |= BAM_PLF_GLF; break;
		case '2': d->format |= BAM_PLF_2ND; break;
		case 'P': d->format |= BAM_PLF_READPOS; break;
		case 'I': d->ido->q_indel = atoi(optarg); break;
		case 'G': d->ido->r_indel = atof(optarg); break;
		case 'S': is_SAM = 1; break;
		case 'R':
			if (strcmp(optarg, "random") == 0) d->format |= BAM_PLF_RANBASE;
			else if (strcmp(optarg, "first") == 0) d->format |= BAM_PLF_1STBASE;
			else if (strcmp(optarg, "all") == 0) d->format |= BAM_PLF_ALLBASE;
			else fprintf(stderr, "[bam_pileup] unrecognized -R\n");
			break;
		default: fprintf(stderr, "Unrecognizd option '-%c'.\n", c); return 1;
		}
	}
	if (fn_list) is_SAM = 1;
	if (optind == argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:  samtools pileup [options] <in.bam>|<in.sam>\n\n");
		fprintf(stderr, "Option: -s        simple (yet incomplete) pileup format\n");
		fprintf(stderr, "        -S        the input is in SAM\n");
		fprintf(stderr, "        -A        use the MAQ model for SNP calling\n");
		fprintf(stderr, "        -2        output the 2nd best call and quality\n");
		fprintf(stderr, "        -i        only show lines/consensus with indels\n");
		fprintf(stderr, "        -m INT    filtering reads with bits in INT [%d]\n", d->mask);
		fprintf(stderr, "        -M INT    cap mapping quality at INT [%d]\n", d->c->cap_mapQ);
        fprintf(stderr, "        -d INT    limit maximum depth for indels [unlimited]\n");
		fprintf(stderr, "        -t FILE   list of reference sequences (force -S)\n");
		fprintf(stderr, "        -l FILE   list of sites at which pileup is output\n");
		fprintf(stderr, "        -f FILE   reference sequence in the FASTA format\n\n");
		fprintf(stderr, "        -c        output the SOAPsnp consensus sequence\n");
		fprintf(stderr, "        -v        print variants only (for -c)\n");
		fprintf(stderr, "        -g        output in the GLFv3 format (suppressing -c/-i/-s)\n");
		fprintf(stderr, "        -T FLOAT  theta in maq consensus calling model (for -c/-g) [%f]\n", d->c->theta);
		fprintf(stderr, "        -N INT    number of haplotypes in the sample (for -c/-g) [%d]\n", d->c->n_hap);
		fprintf(stderr, "        -r FLOAT  prior of a difference between two haplotypes (for -c/-g) [%f]\n", d->c->het_rate);
		fprintf(stderr, "        -G FLOAT  prior of an indel between two haplotypes (for -c/-g) [%f]\n", d->ido->r_indel);
		fprintf(stderr, "        -I INT    phred prob. of an indel in sequencing/prep. (for -c/-g) [%d]\n", d->ido->q_indel);
		fprintf(stderr, "\n");
		free(fn_list); free(fn_fa); free(d);
		return 1;
	}
	if (d->format & (BAM_PLF_RANBASE|BAM_PLF_1STBASE|BAM_PLF_ALLBASE)) d->format |= BAM_PLF_CNS;
	if (fn_fa) d->fai = fai_load(fn_fa);
	if (d->format & (BAM_PLF_CNS|BAM_PLF_GLF)) bam_maqcns_prepare(d->c); // consensus calling
	if (d->format & BAM_PLF_GLF) { // for glf output
		glf3_header_t *h;
		h = glf3_header_init();
		d->fp_glf = bgzf_fdopen(fileno(stdout), "w");
		glf3_header_write(d->fp_glf, h);
		glf3_header_destroy(h);
	}
	if (d->fai == 0 && (d->format & (BAM_PLF_CNS|BAM_PLF_INDEL_ONLY)))
		fprintf(stderr, "[bam_pileup] indels will not be called when -f is absent.\n");
	if (fn_fa && is_SAM && fn_list == 0) fn_list = samfaipath(fn_fa);

	{
		samfile_t *fp;
		fp = is_SAM? samopen(argv[optind], "r", fn_list) : samopen(argv[optind], "rb", 0);
		if (fp == 0 || fp->header == 0) {
			fprintf(stderr, "[bam_pileup] fail to read the header: non-exisiting file or wrong format.\n");
			return 1;
		}
		d->h = fp->header;
		if (fn_pos) d->hash = load_pos(fn_pos, d->h);
		sampileup(fp, d->mask, pileup_func, d);
		samclose(fp); // d->h will be destroyed here
	}

	// free
	if (d->format & BAM_PLF_GLF) bgzf_close(d->fp_glf);
	if (fn_pos) { // free the hash table
		khint_t k;
		for (k = kh_begin(d->hash); k < kh_end(d->hash); ++k)
			if (kh_exist(d->hash, k)) free(kh_val(d->hash, k));
		kh_destroy(64, d->hash);
	}
	free(fn_pos); free(fn_list); free(fn_fa);
	if (d->fai) fai_destroy(d->fai);
	bam_maqcns_destroy(d->c);
	free(d->ido); free(d->ref); free(d);
	return 0;
}

/***********
 * mpileup *
 ***********/

typedef struct {
	char *reg;
	faidx_t *fai;
} mplp_conf_t;

typedef struct {
	bamFile fp;
	bam_iter_t iter;
} mplp_aux_t;

static int mplp_func(void *data, bam1_t *b)
{
	mplp_aux_t *ma = (mplp_aux_t*)data;
	if (ma->iter) return bam_iter_read(ma->fp, ma->iter, b);
	return bam_read1(ma->fp, b);
}

static int mpileup(mplp_conf_t *conf, int n, char **fn)
{
	mplp_aux_t **data;
	int i, tid, pos, *n_plp, beg0 = 0, end0 = 1u<<29, ref_len, ref_tid;
	const bam_pileup1_t **plp;
	bam_mplp_t iter;
	bam_header_t *h = 0;
	char *ref;
	// allocate
	data = calloc(n, sizeof(void*));
	plp = calloc(n, sizeof(void*));
	n_plp = calloc(n, sizeof(int*));
	// read the header and initialize data
	for (i = 0; i < n; ++i) {
		bam_header_t *h_tmp;
		data[i] = calloc(1, sizeof(mplp_aux_t));
		data[i]->fp = bam_open(fn[i], "r");
		h_tmp = bam_header_read(data[i]->fp);
		if (conf->reg) {
			int beg, end;
			bam_index_t *idx;
			idx = bam_index_load(fn[i]);
			if (idx == 0) {
				fprintf(stderr, "[%s] fail to load index for %d-th input.\n", __func__, i+1);
				exit(1);
			}
			if (bam_parse_region(h_tmp, conf->reg, &tid, &beg, &end) < 0) {
				fprintf(stderr, "[%s] malformatted region or wrong seqname for %d-th input.\n", __func__, i+1);
				exit(1);
			}
			if (i == 0) beg0 = beg, end0 = end;
			data[i]->iter = bam_iter_query(idx, tid, beg, end);
			bam_index_destroy(idx);
		}
		if (i == 0) h = h_tmp;
		else {
			// FIXME: to check consistency
			bam_header_destroy(h_tmp);
		}
	}
	// mpileup
	ref_tid = -1; ref = 0;
	iter = bam_mplp_init(n, mplp_func, (void**)data);
	while (bam_mplp_auto(iter, &tid, &pos, n_plp, plp) > 0) {
		if (conf->reg && (pos < beg0 || pos >= end0)) continue; // out of the region requested
		if (tid != ref_tid) {
			free(ref);
			if (conf->fai) ref = fai_fetch(conf->fai, h->target_name[tid], &ref_len);
			ref_tid = tid;
		}
		printf("%s\t%d\t%c", h->target_name[tid], pos + 1, (ref && pos < ref_len)? ref[pos] : 'N');
		for (i = 0; i < n; ++i) {
			int j;
			printf("\t%d\t", n_plp[i]);
			if (n_plp[i] == 0) printf("*\t*");
			else {
				for (j = 0; j < n_plp[i]; ++j)
					pileup_seq(plp[i] + j, pos, ref_len, ref);
				putchar('\t');
				for (j = 0; j < n_plp[i]; ++j) {
					const bam_pileup1_t *p = plp[i] + j;
					int c = bam1_qual(p->b)[p->qpos] + 33;
					if (c > 126) c = 126;
					putchar(c);
				}
			}
		}
		putchar('\n');
	}
	bam_mplp_destroy(iter);
	bam_header_destroy(h);
	for (i = 0; i < n; ++i) {
		bam_close(data[i]->fp);
		if (data[i]->iter) bam_iter_destroy(data[i]->iter);
		free(data[i]);
	}
	free(data); free(plp); free(ref); free(n_plp);
	return 0;
}

int bam_mpileup(int argc, char *argv[])
{
	int c;
	mplp_conf_t mplp;
	memset(&mplp, 0, sizeof(mplp_conf_t));
	while ((c = getopt(argc, argv, "f:r:")) >= 0) {
		switch (c) {
		case 'f':
			mplp.fai = fai_load(optarg);
			if (mplp.fai == 0) return 1;
			break;
		case 'r': mplp.reg = strdup(optarg);
		}
	}
	if (argc == 1) {
		fprintf(stderr, "Usage: samtools mpileup [-r reg] [-f in.fa] in1.bam [in2.bam [...]]\n");
		return 1;
	}
	mpileup(&mplp, argc - optind, argv + optind);
	free(mplp.reg);
	if (mplp.fai) fai_destroy(mplp.fai);
	return 0;
}
