#!/usr/bin/env python3
## ecolityping.py
## Ecoli related analysis
## Runs KMA (k-mer alignment) of reads against a database of genes.


import argparse
import sys
import os


parser = argparse.ArgumentParser(description='ecolityping.py')
parser.add_argument("-i","--sampleid", help='ecoli1')
parser.add_argument("-R1", "--read1", help="ecoli1_R1.fastq.gz")
parser.add_argument("-R2", "--read2", help="ecoli2_R2.fastq.gz")
parser.add_argument("-db", "--db_path", default="db")
parser.add_argument("-k", "--kma_path", default='/usr/bin/')
parser.add_argument("-o", "--output_dir", default='output')
parser.add_argument("-u", "--update", help="--update:yes/no")
args=parser.parse_args()


# Functions
def make_folder_if_not_exists(folder_name):
	"""This function creates output folders if they don't exists."""
	if not os.path.exists(folder_name):
		os.makedirs(folder_name)


# Main
# Update dbs
ECOLIGENESDB = os.path.abspath(os.path.join(args.db_path, "ecoligenes", "ecoligenes"))
if args.update == "yes" or args.update == "True":
	print('# Updating db...')
	# Run KMA
	command=f"{args.kma_path}/kma_index -i {ECOLIGENESDB}.fsa -o {args.db_path}/ecoligenes/ecoligenes"
	print(f"Command {command}") 
	os.system(command)


# Check if any of the optional arguments were provided
optional_args = [args.sampleid, args.read1, args.read2, args.kma_path, args.output_dir]
if all(optional_args):
    print("# Arguments provided...")
else:
    print("# No further task to be done...")
    sys.exit(1)


# Create outputs
SPECOLIFBIDIR=os.path.join(args.output_dir, args.sampleid, "sp_ecoli_fbi") # KMA results
make_folder_if_not_exists(SPECOLIFBIDIR)


# Run KMA
command=f"{args.kma_path}/kma -ipe {args.read1} {args.read2} -matrix -t_db {ECOLIGENESDB} -o {SPECOLIFBIDIR}/colipost"
print(f"Command {command}") 
os.system(command)
