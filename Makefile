# Makefile

.PHONY: test update_db type process summarize

clean:	
	rm -rf test/SRR26510933/

# Combined test target to run all steps
# test: serotyping processing summarizing
test: serotyping md5serotypecheck

# Step 1: Run the typing process
serotyping:
	@BIN_PATH=$$(dirname $$(which python)) && \
	echo "Conda environment bin path: $$BIN_PATH" && \
	python3 ecolityping.py -i SRR26510933 -R1 test/SRR26510933_1_subset.fastq.gz -R2 test/SRR26510933_2_subset.fastq.gz -db db -k "$$BIN_PATH" -o test --update no

# Step 2: Check MD5 checksum
md5serotypecheck:
	cd test && \
	gzip -d SRR26510933/sp_ecoli_fbi/*.gz && \
	md5sum -c Serotypesums.md5
 
# Step 3: Run the post-processing
processing:
	python3 ecolityping.py -i SRR26510933 -R1 test/SRR26510933_1.fastq.gz -R2 test/SRR26510933_2.fastq.gz -db db -k "$$BIN_PATH" -o test --update no

# Step 4: Check MD5 checksum
md5processcheck:
	cd test && \
	md5sum -c Processsums.md5

# Step 5: Run the post-processing
processing:
	python3 ecolityping.py -i SRR26510933 -R1 test/SRR26510933_1.fastq.gz -R2 test/SRR26510933_2.fastq.gz -db db -k "$$BIN_PATH" -o test --update no

# Step 6: Check MD5 checksum
md5processcheck:
	cd test && \
	md5sum -c Summarysums.md5

