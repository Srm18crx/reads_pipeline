# ONT-reads
raw reads analyses

input files are in the fastq_pass folder. The sample used here is 'edna_549_1'
### reads are a mix of 4 different amplicons COI (365), mifish (217 bp), fish_2kb (1974 bp), metazoan_2kb (2290 bp)
NanoPlot, NanoFilt, Seqkit, NGSpeciesID were installed using conda 
## 1. to merge all the fastq.gz files in the fatq_pass folder and gunzip 
cat fastq_pass/*.fastq.gz > all_reads_549_1.fastq.gz

gunzip all_reads_549_1.fastq.gz

## 2. quality check (Nanoplot https://github.com/wdecoster/NanoPlot)
NanoPlot --fastq all_reads_549_1.fastq -o output_raw_reads

## 3. sorting the reads by length in different files (https://github.com/wdecoster/nanofilt)
NanoFilt ../all_reads_549_1.fastq -l 198 --maxlength 321 > mifish_reads_length.fastq

NanoFilt ../all_reads_549_1.fastq -l 345 --maxlength 500 > coi_reads_length.fastq

NanoFilt ../all_reads_549_1.fastq -l 1800 --maxlength 2100 > fish_2_reads_length.fastq

NanoFilt ../all_reads_549_1.fastq -l 1400 --maxlength 2600 > metazoan_reads_length.fastq

### check how many reads are present in each group (Seqkit https://bioinf.shenwei.me/seqkit/download/)
Seqkit stats *_reads_length.fastq

## 4. trimming adapters and primers with cutadapt
mkdir ../cut_primers

cd ../cut_primers
### mifish
cutadapt -g GTCGGTAAAACTCGTGCCAGC...CAAACTGGGATTAGATACCCCACTATG -o mifish_549_trimmed.fastq ../raw_reads/mifish_reads_length.fastq -e 0.2 -m 150 -M 200 --discard-untrimmed

cutadapt -g CATAGTGGGGTATCTAATCCCAGTTTG...GCTGGCACGAGTTTTACCGAC -o mifish_549_trimmed_reverse.fastq ../raw_reads/mifish_reads_length.fastq -e 0.2 -m 150 -M 200 --discard-untrimmed

mkdir final_reads

cat mifish_*.fastq > final_reads/reads_mifish.fastq

### coi
cutadapt -g GGWACWRGWTGRACWNTNTAYCCYCC...TGRTTYTTYGGNCAYCCNGARGTNTA -o coi_549_trimmed.fastq ../raw_reads/coi_reads_length.fastq -e 0.2 -m 345 -M 385 --discard-untrimmed

cutadapt -g TANACYTCNGGRTGNCCRAARAAYCA...GGRGGRTANANWGTYCAWCYWGTWCC -o coi_549_trimmed_reverse.fastq ../raw_reads/coi_reads_length.fastq -e 0.2 -m 345 -M 385 --discard-untrimmed

cat coi_*.fastq > final_reads/reads_coi.fastq

### metazoan
cutadapt -g GTGCCAGCHNHHGCGGTYA...ACRTGAKYTGAGTTCARAYCGG -o metazoan_549_trimmed.fastq ../raw_reads/metazoan_reads_length.fastq -e 0.2 -m 1400 -M 2500 --discard-untrimmed

cutadapt -g CCGRTYTGAACTCARMTCAYGT...TRACCGCDDNDGCTGGCAC -o metazoan_549_trimmed_reverse.fastq ../raw_reads/metazoan_reads_length.fastq -e 0.2 -m 1400 -M 2500 --discard-untrimmed

cat metazoan_*.fastq > final_reads/reads_metazoan.fastq

### fish2 kb
cutadapt -g GGATTAGATACCCYACTATGC...CTAGGGATAACAGCGCAATC -o fish_549_trimmed.fastq ../raw_reads/fish_2_reads_length.fastq -e 0.2 -m 1900 -M 2083 --discard-untrimmed

cutadapt -g GATTGCGCTGTTATCCCTAG...GCATAGTRGGGTATCTAATCC -o fish_549_trimmed_reverse.fastq ../raw_reads/fish_2_reads_length.fastq -e 0.2 -m 1800 -M 2083 --discard-untrimmed

cat fish_*.fastq > final_reads/reads_fish.fastq

## 5. building consensus sequences 
NGSpeciesID --ont --fastq ../reads_metazoan.fastq  --consensus --medaka --abundance_ratio 0.00005 --q 10 --rc_identity_threshold 0.97 --max_seqs_for_consensus 200 --outfolder ./output_metazoan_00005_97_200 --t 48

## 6. blast consensus sequences
blastn -query metazoan_consensus.fasta -out output_blast_metazoan.txt -db ~/Desktop/Sara/NCBI_nt_db/nt -outfmt '6 qseqid sseqid pident evalue length mismatch gapopen gaps qstart qend staxids sscinames scomnames sskingdoms stitle' -perc_identity 90 -max_target_seqs 1 -num_threads 24
