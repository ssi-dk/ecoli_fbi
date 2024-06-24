#!/usr/bin/env python3
## postecolityping.py
## Ecoli related analysis
## Parses ecolityping.sh results to generate an assessment report.


import argparse
import os
import time
import datetime
import json


parser = argparse.ArgumentParser(description='postecolityping.py')
parser.add_argument("-i","--sampleid", help='ecoli1')
parser.add_argument("-d", "--indir")
parser.add_argument("-stbit", help="ST string", default='ST:NA,NA')
args=parser.parse_args()


# Functions
def get_rundir(indir) -> str:
	"""This function gets the rundir from a dir path."""
	rundir = os.path.basename(indir)
	return rundir


def get_folder_creation_date(folder_path):
    creation_time = os.path.getctime(folder_path)
    creation_date = datetime.datetime.fromtimestamp(creation_time)
    formatted_date = str(creation_date.strftime('%y%m%d'))
    return formatted_date


def convert_number_to_date(number):
	"""This function converts a number to a date from 221117 to a 2022-11-17."""
	year=int(str(number)[:2])
	month=int(str(number)[2:4])
	day=int(str(number)[4:6])
	# Create a date string in the "YYYY-MM-DD" format
	date_string=f"20{year:02d}-{month:02d}-{day:02d}"
	return date_string


def get_wgs_date_and_number(rundir):
	"""This function separates the rundir into date and experiment_name/wgsnumber."""
	if "_" in rundir:
		#231006_NB551234_0051_N_WGS_743_AHNLHHAFX5
		rundir_list=rundir.split("_") 	# ['231006', 'NB551234', '0051', 'N', 'WGS', '743', 'AHNLHHAFX5']
		wgsdate = convert_number_to_date(rundir_list[0]) # 231006 -> 2023-10-06
		wgsnumber="_".join(rundir_list[3:6]) # N_WGS_743
	else:
		wgsnumber = rundir
		date = get_folder_creation_date(rundir)
		wgsdate = convert_number_to_date(date)
	return wgsdate, wgsnumber


def print_header_to_output(OUTFILE):
	header="isolate\twzx\twzy\tfliC\tOH\tstx\teae\tehx\tother\twgsrun\twgsdate\tST\tSTgenes\tverbose\n"
	outfile=open(OUTFILE, 'w')
	print(header, end='', file=outfile)
	return header


def write_to_json(csv_data, json_outfile):
	"""Write to json file
	"""
	with open(json_outfile, "w") as jsonout:
		json.dump(csv_data, jsonout)


# Main
# KMA results processing
SPECOLIFBIDIR=os.path.join(args.indir, args.sampleid, "sp_ecoli_fbi") # KMA results
OUTFILE=os.path.join(args.indir, args.sampleid, f"{args.sampleid}.tsv")
rundir = get_rundir(args.indir)
wgsdate, wgsnumber = get_wgs_date_and_number(args.indir)
ST_list=args.stbit.split(",") # ['ST:11', 'adk:12', 'fumC:12', 'gyrB:8', 'icd:12', 'mdh:15', 'purA:2', 'recA:2']
ST=ST_list.pop(0).split(":")[1] # 11
STgenes=",".join(ST_list) # adk:12,fumC:12,gyrB:8,icd:12,mdh:15,purA:2,recA:2
OUTFILE=os.path.join(args.indir, args.sampleid, f"{args.sampleid}.tsv")
header = print_header_to_output(OUTFILE)


# Settings
fliflag="false"
locations={"OH": 4, "stx":5, "wzx":1, "wzy":2, "wzt":1, "wzm":2, "fliC":3, "fli":3, "eae": 6, "ehxA":7, "wgsrun":8, "wgsdate":9, "ST": 10, "STgenes":11, "other":12}
thresholds={"stx":[98,98], "wzx":[98,98], "wzy":[98,98], "wzt":[98,98], "wzm":[98,98], "fliC":[90,90], "fli":[90,90], "eae": [95,95], "ehxA":[95,95], "other":[98,98]}
stxbadthreshold=[30,90]


# Processing of .res file
hitdict={}
hitdict[args.sampleid]={}
res_file=open(os.path.join(SPECOLIFBIDIR, "colipost.res"), 'r')
for line in res_file:
	if line.startswith("#"):
		continue
	#line[0]=the specific hit, line[9]=pc_idt for the hit to the ref, line[18]=the mutation of the specific line
	line=line.split("\t") 
	template=line[0]
	template_cov=line[5]
	query_ident=line[6]
	gene_list=template.split("__")
	gene_name=gene_list[1]
	serotype_or_virulence=gene_list[2]
	if gene_name.startswith("stx"):
		gene_name="stx"
	if not gene_name in hitdict[args.sampleid]:
		hitdict[args.sampleid][gene_name]=[]#{"lencov":lencov, "cov":line[5], "SNP":0, "DEL":0, "INS":0, ".":0}
	hit_results_list=[serotype_or_virulence, template_cov, query_ident]
	hitdict[args.sampleid][gene_name].append(hit_results_list)
#print('HITDICT',hitdict)


# Processing of hits
results_dict={}
for args.sampleid in hitdict:
	foundgoodstx="false"
	results_dict[args.sampleid]=[[],[],[],[],[],[],[],[],[],[],[],[],[],[]]# Remember to increase this if "locations" is made longer  [[], [], [],[],[],[],[],[],[],[],[],,[],[][]]
	Opick="-"	
	# First the mapped hits are parsed and types are decided upon based on a logic from KNJ and SSC
	for gene_name in hitdict[args.sampleid]:					
		for hit_results_list in hitdict[args.sampleid][gene_name]:
			serotype_or_virulence=hit_results_list[0]
			template_cov=hit_results_list[1]
			query_ident=hit_results_list[2]
			if gene_name.startswith("fli") and not "fliC" in gene_name and fliflag=="false":
				results_dict[args.sampleid][locations["fliC"]]=[]
				fliflag="true"
			if not gene_name in thresholds:
				gene_name="other"
			if fliflag=="true":
				gene_name="fli"
			if "stx" in gene_name:
				serotype_or_virulence=serotype_or_virulence.lower()
			if "100.00" in template_cov and "100.00" in query_ident:
				if "wzt" in gene_name or "wzm" in gene_name or fliflag=="true":
					results_dict[args.sampleid][-1].append(":".join(["_".join([gene_name, serotype_or_virulence]), template_cov.lstrip(), query_ident.lstrip()]))
				if gene_name.startswith("e"):
					serotype_or_virulence="positive"
				if gene_name.startswith("stx"):
					foundgoodstx="true"
			elif float(template_cov)> thresholds[gene_name][0] and float(query_ident)>thresholds[gene_name][1]:
				if gene_name.startswith("wz") or fliflag=="true":
					results_dict[args.sampleid][-1].append(":".join(["_".join([gene_name, serotype_or_virulence]), template_cov.lstrip(), query_ident.lstrip()]))
				else:
					results_dict[args.sampleid][-1].append(":".join([serotype_or_virulence, template_cov.lstrip(), query_ident.lstrip()]))
				if gene_name.startswith("e"):
					serotype_or_virulence="positive"
				if gene_name.startswith("stx"):
					foundgoodstx="true"
			elif gene_name.startswith("stx") and float(template_cov)>stxbadthreshold[0] and float(query_ident)>stxbadthreshold[1]:
				results_dict[args.sampleid][-1].append(":".join([serotype_or_virulence, template_cov.lstrip(), query_ident.lstrip()]))
				serotype_or_virulence=serotype_or_virulence[:4]
			else:
				if gene_name.startswith("wz"):
					results_dict[args.sampleid][-1].append(":".join(["_".join([gene_name, serotype_or_virulence]), template_cov.lstrip(), query_ident.lstrip()]))
				else:
					results_dict[args.sampleid][-1].append(":".join([serotype_or_virulence, template_cov.lstrip(), query_ident.lstrip()]))
				continue
			if "fliC" in gene_name and fliflag=="true":
				continue
			if not serotype_or_virulence in results_dict[args.sampleid][locations[gene_name]]:
				results_dict[args.sampleid][locations[gene_name]].append(serotype_or_virulence)
	#print(foundgoodstx)
	if foundgoodstx=="true":
		for stxcase in results_dict[args.sampleid][locations["stx"]]:
			if len(stxcase)==4:
				results_dict[args.sampleid][locations["stx"]].pop(results_dict[args.sampleid][locations["stx"]].index(stxcase))
	results_dict[args.sampleid][locations["stx"]]="; ".join(results_dict[args.sampleid][locations["stx"]]).replace("-", "")
	for gene in results_dict[args.sampleid]:
		if type(gene)==list:
			results_dict[args.sampleid][results_dict[args.sampleid].index(gene)]=":".join(gene)
	if len(results_dict[args.sampleid][locations["wzx"]]) > 1 and not "*" in results_dict[args.sampleid][locations["wzx"]]:
		Opick=results_dict[args.sampleid][locations["wzx"]]
	if len(results_dict[args.sampleid][locations["wzy"]]) > 1 and not "*" in results_dict[args.sampleid][locations["wzy"]]:
		Opick=results_dict[args.sampleid][locations["wzy"]]
	
	if len(results_dict[args.sampleid][locations["fliC"]]) < 2:
		results_dict[args.sampleid][locations["fliC"]]="-"
	results_dict[args.sampleid][locations["OH"]]=":".join([Opick, results_dict[args.sampleid][locations["fliC"]]])

	# Here the results_dict is curated with dashes
	for element in results_dict[args.sampleid]:
		if len(element)<2:
			results_dict[args.sampleid][results_dict[args.sampleid].index(element)]="-"
	results_dict[args.sampleid].pop(0)
	results_dict[args.sampleid][locations["wgsrun"]]=wgsnumber
	#print(wgsnumber)
	results_dict[args.sampleid][locations["wgsdate"]]=wgsdate
#print(results_dict)	


# Here the results_dict is curated with dashes
for args.sampleid in results_dict.keys():
	outfile=open(OUTFILE, 'w')
	results_dict[args.sampleid][locations["ST"]]=ST # as defined on line 34 KRKI 02-05-2019
	results_dict[args.sampleid][locations["STgenes"]]=STgenes # as defined on line 35 KRKI 02-05-2019
	lineelements=[]
	for element in results_dict[args.sampleid]:
		if type(element) is dict:
			for gene in ["O", "H", "stx", "eae", "ehx"]:
				if gene in element:
					lineelements.append(";".join(element[gene]))
		elif len(element)<1:
			#print('element', element)
			lineelements.append("-")
		else:
			lineelements.append(element)

	
	# Print results to outfile
	print("".join([header, "\t".join([args.sampleid, "\t".join(lineelements)])]), file=outfile, end='')


# Turn the data into a dictionary
keys = header.strip().split('\t')
csv_data = {key: None for key in keys}
lineelements.insert(0, args.sampleid)
csv_data = dict(zip(keys, lineelements))
#print(csv_data)
json_outfile=os.path.join(os.path.join(args.indir, args.sampleid, f"{args.sampleid}.json"))
write_to_json(csv_data, json_outfile)
