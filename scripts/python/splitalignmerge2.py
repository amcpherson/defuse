import cmdrun
import optparse
import csv
import os
import sys

supported_aligners = ["bowtie","mrsfast"]
supported_align_types = ["paired","single"]

if __name__ == '__main__':

	usage = "usage: \%prog [options] reads_end_1_fastq reads_end_2_fastq [reference_fasta output_bam]"

	parser = optparse.OptionParser(usage)
	parser.add_option("-w", "--working_dir", type="string", dest="working_dir", help="Working directory for temp files", default="./")
	parser.add_option("-n", "--num_job_reads", type="string", dest="num_job_reads", help="Number of reads per split", default="1000000")
	parser.add_option("-a", "--aligner", type="choice", dest="aligner", help="Alignment tool", default=supported_aligners[0], choices=supported_aligners)
	parser.add_option("-t", "--align_type", type="choice", dest="align_type", help="Alignment type, either paired or single", default=supported_align_types[0], choices=supported_align_types)
	parser.add_option("-b", "--aligner_bin", type="string", dest="aligner_bin", help="Path to alignment tool binary", default="")
	parser.add_option("-s", "--samtools_bin", type="string", dest="samtools_bin", help="Path to samtools binary", default="samtools")
	parser.add_option("-f", "--fragment_length", type="string", dest="fragment_length", help="Maximum fragment length", default="500")
	parser.add_option("-m", "--max_alignments", type="string", dest="max_alignments", help="Maximum alignments per read", default="1")
	parser.add_option("-e", "--edit_distance", type="string", dest="edit_distance", help="Maximum edit distance for a valid alignment", default="2")
	parser.add_option("-r", "--remove_split", type="string", dest="remove_split", help="Remove split files when finished", default=False)

	(options, args) = parser.parse_args()

	# Aligner binary filename defaults to name of aligner
	if options.aligner_bin == "":
		options.aligner_bin = options.aligner

	# Two fastq files required
	if len(args) < 2:
		parser.error("incorrect number of arguments")
	reads_end_1_fastq = os.path.abspath(args[0])
	reads_end_2_fastq = os.path.abspath(args[1])
	
	splitalignmerge(*args, **options)

def splitalignmerge(reads_end_1_fastq, reads_end_2_fastq, reference_fasta="", output_bam="",
                    working_dir="./", num_job_reads="1000000", aligner=supported_aligners[0],
                    align_type=supported_align_types[0], aligner_bin="", samtools_bin="samtools",
                    fragment_length="500", max_alignments="1", edit_distance="2", split_only=True,
                    remove_split=False):

	reads_end_1_fastq = os.path.abspath(reads_end_1_fastq)
	reads_end_2_fastq = os.path.abspath(reads_end_2_fastq)

	# Only split if no reference or output specified
	split_only = True
	if len(args) == 4:
		reference_fasta = os.path.abspath(reference_fasta)
		output_bam = os.path.abspath(output_bam)
		split_only = False
	elif len(args) != 2:
		parser.error("incorrect number of arguments")
	
	# Fasta index is required for samtools
	reference_fasta_index = reference_fasta + ".fai";
	if not os.path.exists(reference_fasta_index):
		sys.stderr.write("Error: Required file " + reference_fasta_index + " does not exist\n")
		sys.exit(1)
	
	# Ensure everything is absolute paths
	options.working_dir = os.path.abspath(options.working_dir)
	
	# Use the reads basename and reference basename for uniquely naming files
	reads_name = os.path.basename(reads_end_1_fastq)
	reference_name = os.path.basename(reference_fasta)
	
	name = reads_name + "." + reference_name
	
	cmdrun = cmdrun.cmdrun(name, options.working_dir, "sge")
	
	split_prefix = options.working_dir + "/" + reads_name + "." + options.num_job_reads
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
	
	def align(fastq_end_1, fastq_end_2, bam_prefix, bam_filenames):
	
		align_command = []
	
		if options.aligner == "bowtie":
			align_command.extend([options.aligner_bin, "--sam-nosq", "-S", "--mm", "-t", "-k", options.max_alignments, "-m", options.max_alignments, reference_fasta])
			if options.align_type == "paired":
				align_command.extend(["-X", options.fragment_length, "-1", "%(fastq1)s", "-2", "%(fastq2)s"])
			elif options.align_type == "single":
				align_command.extend(["--min", "0", "--max", options.fragment_length, "%(fastq)s"])
	
		elif options.aligner == "mrsfast":
			align_command.extend([options.aligner_bin, "-e", options.edit_distance, "-n", options.max_alignments, "--search", reference_fasta])
			if options.align_type == "paired":
				align_command.extend(["--min", "0", "--max", options.fragment_length, "-seq1", "%(fastq1)s", "-seq2", "%(fastq2)s"])
			elif options.align_type == "single":
				align_command.extend(["-seq", "%(fastq)s"])
			
		bam_sort_prefix = bam_prefix + ".sort"
		
		align_command.extend(["|", options.samtools_bin, "view", "-bt", reference_fasta_index, "-", "|", options.samtools_bin, "sort", "-o", "-", bam_sort_prefix, ">", "%(bam)s"])
	
		if options.align_type == "paired":
			bam_filename = bam_prefix + ".bam"
			cmdrun.padd(align_command, {'fastq1': fastq_end_1, 'fastq2': fastq_end_2}, {'bam': bam_filename})
			bam_filenames.append(bam_filename)
	
		elif options.align_type == "single":
			bam_filename = bam_prefix + ".1.bam"
			bam_filename = bam_prefix + ".2.bam"
			cmdrun.padd(align_command, {'fastq': fastq_end_1}, {'bam': bam_filename})
			cmdrun.padd(align_command, {'fastq': fastq_end_2}, {'bam': bam_filename})
			bam_filenames.append(bam_filename_1)
			bam_filenames.append(bam_filename_2)
	
	
	def remove_files(filenames):
		for filename in filenames:
			if os.path.exists(filename):
				os.remove(filename)
	
	if split_only:
		if options.remove_split:
			remove_files(split_filenames)
			os.remove(split_catalog_filename)
		sys.exit(0)
	
	# Create bam files
	bam_filenames = []
	for split_filename_pair in split_filename_pairs:
		reads_end_1_split_fastq = split_filename_pair[0]
		reads_end_2_split_fastq = split_filename_pair[1]
		split_name = os.path.basename(reads_end_1_split_fastq)
		bam_prefix = options.working_dir + "/" + split_name + "." + reference_name
		align(reads_end_1_split_fastq, reads_end_2_split_fastq, bam_prefix, bam_filenames)
	
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
			intermediate_bam_filename = options.working_dir + "/" + name + "." + reference_name + ".intermediate.%d.bam" % len(intermediate_bam_filenames)
	
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
	
