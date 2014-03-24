import cmdrun
import optparse
import csv
import os
import sys
import hashlib

usage = """usage: %prog [options] input output command

command must include %(input)s and %(output)s""";

parser = optparse.OptionParser(usage)
parser.add_option("-s", "--split_size", type="int", dest="split_size", help="Split size in lines", default="1000")
parser.add_option("-w", "--working_dir", type="string", dest="working_dir", help="Working directory", default="./")
parser.add_option("-t", "--submit_type", type="string", dest="submit_type", help="Type of job submission", default="direct")
parser.add_option("-m", "--job_mem", type="string", dest="job_mem", help="Memory per job in Gigabytes", default="2")

(options, args) = parser.parse_args()

if len(args) < 3:
	parser.error("incorrect number of arguments")

input = os.path.abspath(args[0])
output = os.path.abspath(args[1])
command = args[2]

options.working_dir = os.path.abspath(options.working_dir)

m = hashlib.md5()
m.update(input.encode("utf-8"))
m.update(output.encode("utf-8"))
m.update(command.encode("utf-8"))
job_name = "mapmerge." + m.hexdigest()

print(job_name)

log_filename = options.working_dir + "/" + job_name + ".log"

cmdrun = cmdrun.cmdrun(job_name, options.working_dir, log_filename, job_type=options.submit_type, job_mem=options.job_mem)

def remove_files(filenames):
	for filename in filenames:
		if os.path.exists(filename):
			os.remove(filename)

split_input_filenames = []
split_input_file = 0
split_input_index = 0
split_input_lines = options.split_size

input_file = open(input, 'r')
for line in input_file:
	if (split_input_lines == options.split_size):
		if (split_input_file):
			split_input_file.close()
		split_input_filename = options.working_dir + "/" + job_name + ".%d" % split_input_index
		split_input_filenames.append(split_input_filename)
		split_input_file = open(split_input_filename, 'w')
		split_input_index += 1
		split_input_lines = 0
	split_input_file.write(line)
	split_input_lines += 1

if (split_input_file):
	split_input_file.close()

split_output_filenames = []
for split_input_filename in split_input_filenames:
	split_output_filename = split_input_filename + ".out"
	split_output_filenames.append(split_output_filename)
	cmdrun.padd([command], {'input': split_input_filename}, {'output': split_output_filename});
cmdrun.prun()

merge_command = ["cat"]
merge_inputs = {}
for merge_index, split_output_filename in enumerate(split_output_filenames):
	merge_input_id = "merge%d" % merge_index
	merge_command.append("%(" + merge_input_id + ")s")
	merge_inputs[merge_input_id] = split_output_filename
merge_command.extend([">", "%(output)s"])
cmdrun.run(merge_command, merge_inputs, {'output': output})

remove_files(split_input_filenames)
remove_files(split_output_filenames)
os.remove(log_filename)

