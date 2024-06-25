#!/usr/bin/env python3
## Ecoli related analysis
## Summarizes the contents of various .tsv files by keeping 1 header and concatenating the contents.


import argparse
import os
import subprocess


parser = argparse.ArgumentParser(description='qcecolisummary.py')
parser.add_argument("-i", "--input_dir", default='input')
parser.add_argument("-o", "--output_dir", default='output')
args = parser.parse_args()


# Functions
def make_folder_if_not_exists(folder_name):
	"""This function creates output folders if they don't exists."""
	if not os.path.exists(folder_name):
		os.makedirs(folder_name)


def find_dirs(path):
	"""This function finds dirs in a path,
	excluding files."""
	dirs = [os.path.join(path, dir) for dir in os.listdir(path) if os.path.isdir(os.path.join(path, dir))]
	return dirs


def find_files(path, extension=None):
    """This function finds files in a path with optional extension filtering,
    excluding directories."""
    if extension:
        files = [os.path.join(path, file) for file in os.listdir(path)
                 if os.path.isfile(os.path.join(path, file)) and file.endswith(extension)]
    else:
        files = [os.path.join(path, file) for file in os.listdir(path)
                 if os.path.isfile(os.path.join(path, file))]
    return files


# Main
# Create outputs
make_folder_if_not_exists(args.output_dir)


sampledirs = find_dirs(args.input_dir)
for sampledir in sampledirs:
	outfile_path = f"{os.path.join(args.output_dir, os.path.basename(args.input_dir))}.tsv"
	if os.path.exists(outfile_path):
		outfile=open(outfile_path,"a")
		header="true"
		line_ending='\n'
	else:
		outfile=open(outfile_path,"w")
		header="false"
	files = find_files(sampledir, '.tsv')
	for file in files:
		print(file)
		if header=="false":
			proc=subprocess.Popen("".join(["head -1 ", file]), stdout=subprocess.PIPE, shell=True)
			(out, err)=proc.communicate()
			outfile.write(out.decode())
			header="true"
		proc=subprocess.Popen("".join(["tail -n +2 ", file]), stdout=subprocess.PIPE, shell=True)
		(out, err)=proc.communicate()
		outfile.write(out.decode() + '\n')
