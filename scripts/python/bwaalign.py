import cmdrun
import optparse
import csv
import os
import sys

usage = "usage: %prog [options] reads_end_1_fastq reads_end_2_fastq reference_fasta output_bam"
parser = optparse.OptionParser(usage)
parser.add_option("-n", "--num_job_reads", type="string", dest="num_job_reads", help="Number of reads per job", default="0")
parser.add_option("-b", "--bwa_bin", type="string", dest="bwa_bin", help="Path to bwa binary", default="bwa")
parser.add_option("-s", "--samtools_bin", type="string", dest="samtools_bin", help="Path to samtools binary", default="samtools")
parser.add_option("-f", "--fragment_length", type="string", dest="fragment_length", help="Maximum fragment length", default="500")
parser.add_option("-r", "--remove_split", type="string", dest="remove_split", help="Remove split files when finished", default=False)
parser.add_option("-w", "--working_dir", type="string", dest="working_dir", help="Working directory for temp files", default="./")
parser.add_option("-t", "--submit_type", type="string", dest="submit_type", help="Type of job submission", default="direct")
parser.add_option("-m", "--max_parallel", type="int", dest="max_parallel", help="Maximum number of parallel jobs", default=200)

(options, args) = parser.parse_args()

# Two fastq files required
if len(args) < 4:
	parser.error("incorrect number of arguments")

reads_end_1_fastq = os.path.abspath(args[0])
reads_end_2_fastq = os.path.abspath(args[1])
reference_fasta = os.path.abspath(args[2])
output_bam = os.path.abspath(args[3])

# Ensure working directory is absolute paths
options.working_dir = os.path.abspath(options.working_dir)

job_name = os.path.basename(reads_end_1_fastq) + "." + os.path.basename(reference_fasta)
log_filename = options.working_dir + "/" + job_name + ".log"
cmdrun = cmdrun.cmdrun(job_name, options.working_dir, log_filename, options.submit_type, options.max_parallel)

def remove_files(filenames):
	for filename in filenames:
		if os.path.exists(filename):
			os.remove(filename)

if options.num_job_reads == "0":
	
	reads_end_1_sai = options.working_dir + "/" + os.path.basename(reads_end_1_fastq) + "." + os.path.basename(reference_fasta) + ".sai"
	reads_end_2_sai = options.working_dir + "/" + os.path.basename(reads_end_2_fastq) + "." + os.path.basename(reference_fasta) + ".sai"
	bam_sort_prefix = output_bam + ".sort"
	
	cmdrun.run([options.bwa_bin, "aln", reference_fasta, "%(fastq)s", ">", "%(sai)s"],
			   {'fastq': reads_end_1_fastq}, {'sai': reads_end_1_sai})
	cmdrun.run([options.bwa_bin, "aln", reference_fasta, "%(fastq)s", ">", "%(sai)s"],
			   {'fastq': reads_end_2_fastq}, {'sai': reads_end_2_sai})
	
	cmdrun.run([options.bwa_bin, "sampe", "-a", options.fragment_length, reference_fasta, "%(sai1)s", "%(sai2)s", "%(fastq1)s", "%(fastq2)s", "|", options.samtools_bin, "view", "-bt", reference_fasta + ".fai", "-", "|", options.samtools_bin, "sort", "-o", "-", bam_sort_prefix, ">", "%(bam)s"],
			   {'sai1': reads_end_1_sai, 'sai2': reads_end_2_sai, 'fastq1': reads_end_1_fastq, 'fastq2': reads_end_2_fastq}, {'bam': output_bam})
			   
	os.remove(reads_end_1_sai)
	os.remove(reads_end_2_sai)
	
else:
	
	split_prefix = options.working_dir + "/" + os.path.basename(reads_end_1_fastq) + "." + options.num_job_reads
	split_catalog_filename = split_prefix + ".split.catalog"
	
	split_fastq_script = os.path.abspath(sys.path[0]) + "/split_fastq.pl"
	
	cmdrun.run([split_fastq_script, "%(fastq1)s", "%(fastq2)s", options.num_job_reads, split_prefix, "> %(catalog)s"],
			   {'fastq1': reads_end_1_fastq, 'fastq2': reads_end_2_fastq}, {'catalog': split_catalog_filename});
	
	split_filenames = []
	split_filename_pairs = []
	for split_fastq_info in csv.reader(open(split_catalog_filename, 'r'), delimiter='\t'):
		split_filenames.append(split_fastq_info[0])
		split_filenames.append(split_fastq_info[1])
		split_filename_pairs.append([split_fastq_info[0], split_fastq_info[1]])

	# Run alignment jobs
	bam_filenames = []
	for split_filename_pair in split_filename_pairs:
		reads_end_1_split_fastq = split_filename_pair[0]
		reads_end_2_split_fastq = split_filename_pair[1]
		output_split_bam = options.working_dir + "/" + os.path.basename(reads_end_1_split_fastq) + "." + os.path.basename(reference_fasta) + ".bam"
		
		cmdrun.padd([sys.executable, os.path.abspath(sys.argv[0]), "%(fastq1)s", "%(fastq2)s", reference_fasta, "%(outbam)s", "-b", options.bwa_bin, "-s", options.samtools_bin, "-f", options.fragment_length, "-w", options.working_dir],
					{'fastq1': reads_end_1_split_fastq, 'fastq2': reads_end_2_split_fastq}, {'outbam': output_split_bam});
					
		bam_filenames.append(output_split_bam)
	
	# Run the alignment commands
	cmdrun.prun()
	
	# Maximum number of files to merge at once
	merge_max = 100
	
	# Merge merge_max at a time
	current_merge_bam_filenames = []
	intermediate_bam_filenames = []
	for bam_filename in bam_filenames:
	
		# Maintain a list of the current bam filenames
		current_merge_bam_filenames.append(bam_filename)
		
		# Merge current list if we have reached the max or the last filename
		if len(current_merge_bam_filenames) == merge_max or bam_filename == bam_filenames[-1]:
			
			# Create intermediate bam filename
			intermediate_bam_filename = options.working_dir + "/" + job_name + ".intermediate.%d.bam" % len(intermediate_bam_filenames)
	
			# Merge bam files or copy if theres only one
			if len(current_merge_bam_filenames) == 1:
				os.rename(current_merge_bam_filenames[0], intermediate_bam_filename)
			else:
				merge_command = [options.samtools_bin, "merge", "%(intermediate)s"]
				merge_inputs = {}
				for merge_index, current_merge_bam_filename in enumerate(current_merge_bam_filenames):
					merge_input_id = "merge%d" % merge_index
					merge_command.append("%(" + merge_input_id + ")s")
					merge_inputs[merge_input_id] = current_merge_bam_filename
				cmdrun.padd(merge_command, merge_inputs, {'intermediate': intermediate_bam_filename})
			
			# Add intermediate filename to list and clear merge filename list
			intermediate_bam_filenames.append(intermediate_bam_filename)
			current_merge_bam_filenames = []
			
	# Run the merge commands
	cmdrun.prun()
	
	# Merge the intermediate bam files or copy if theres only one
	if len(intermediate_bam_filenames) == 1:
		os.rename(intermediate_bam_filenames[0], output_bam)
	else:
		merge_command = [options.samtools_bin, "merge", "%(output)s"]
		merge_inputs = {}
		for merge_index, intermediate_bam_filename in enumerate(intermediate_bam_filenames):
			merge_input_id = "merge%d" % merge_index
			merge_command.append("%(" + merge_input_id + ")s")
			merge_inputs[merge_input_id] = intermediate_bam_filename
		cmdrun.run(merge_command, merge_inputs, {'output': output_bam})
	
	remove_files(bam_filenames)
	remove_files(intermediate_bam_filenames)
	
	if options.remove_split:
		remove_files(split_filenames)
		os.remove(split_catalog_filename)
		

