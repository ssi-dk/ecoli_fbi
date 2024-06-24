# ecoli_fbi
Determines the serotype and virulence in *Escherichia coli* isolates using NGS data (reads) through kmer alignment.

## Quick start
```bash
# Type
python3 ecolityping.py -h
# Process
python3 postecolityping.py -h
# Summarize
python3 qc_ecoli_summary.py -h
```

## Installation
### Source
```bash
# Clone this repo
git clone https://github.com/ssi-dk/ecoli_fbi.git
# Create an environment with the required tools with conda
conda create --name ecoli_pipeline kma python=3.11
# Activate the environment
conda activate ecoli pipeline
# Install pip requirements
pip install openpyxl pandas envyaml pandas requests python-dotenv
```

## Usage
### Example
```bash
# Download data into the test folder
mkdir -p test
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/033/SRR26510933/SRR26510933_1.fastq.gz -P test
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/033/SRR26510933/SRR26510933_2.fastq.gz -P test

```

### Pipeline
```bash
# Build the db
python3 ecolityping.py -db db -update yes
# Type
python3 ecolityping.py -i SRR26510933 -R1 test/SRR26510933_1.fastq.gz -R2 test/SRR26510933_2.fastq.gz -db db -k /usr/bin -o test --update no
# Process
python3 postecolityping.py -i SRR26510933 -d test -stbit "STNA;NA:NA"
# Summarize
python3 qc_ecoli_summary.py -i test -o test
```

## Output
### .tsv
```
isolate	wzx	wzy	fliC	OH	stx	eae	ehx	other	wgsrun	wgsdate	ST	STgenes	verbose
SRR26510933	-	-	H7	-:H7	stx1a; stx2c	positive	positive	-	240624	2024-06-24	NA	-	wzx_O157:96.53:100.00:wzy_O157:92.74:100.00
```
### .json
```
{"isolate": "SRR26510933", "wzx": "-", "wzy": "-", "fliC": "H7", "OH": "-:H7", "stx": "stx1a; stx2c", "eae": "positive", "ehx": "positive", "other": "-", "wgsrun": "test", "wgsdate": "2024-06-24", "ST": "NA", "STgenes": "-", "verbose": "wzx_O157:96.53:100.00:wzy_O157:92.74:100.00"}
```

## Updating the db
```bash
# Build the db
python3 ecolityping.py -db db --update yes
```
