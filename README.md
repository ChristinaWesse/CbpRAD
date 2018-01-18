# *Capsella bursa-pastoris* RAD-Seq

Last modified: January 8th 2018, Christina Wesse

## Primary Analysis
During 6 sessions from Nov 2014 to Dec 2016 tissue from Shepherd's Purse *Capsella bursa-pastoris* (Cbp in the following) was DNA extracted with CTAB and chloroform. Library prep for single end RAD sequencing was performed with 384 individual adaptor barcodes each multiplex and KpnI as cutter. Sequencing was performed in the Genome Center of the Max Planck Institute for Developmental Biology in Tübingen, Germany. DNA extracts and library backup are stored in the Weigel Lab. 

### Sequencing files
Each sequencing lane is a multiplex of up to 384 individually barcoded samples. Single-end ("R1").
Sequencing files are on the biocluster of the UCR (biocluster.ucr.edu) in
`koeniglab/shared/OUTSIDE_SEQUENCE_DATA/CHRISTINA_CBP_RAD`

For better overview, each lane will be named referring to their sample composition:

#### "Nov14":
`/illumina_SN7001143flowcellB_SampleIdNA_RunId0213_LaneId5/lane5_NoIndex_L005_R1_001.fastq.gz`
File size: 17,1 GB
Total number of reads: 176478793
Barcodes used: 384

#### "Am1"* ("America1"):
`/illumina_ST-J00101flowcellA_SampleIdNA_RunId0008_LaneId7/Undetermined_S0_L007_R1_001.fastq.gz`
File size: 19,1 GB
Total number of reads: 323376947
Barcodes used: 232

#### "Am2"* ("America2"):
`/illumina_ST-J00101flowcellA_SampleIdNA_RunId0008_LaneId8/Undetermined_S0_L008_R1_001.fastq.gz`
File size: 18,3 GB
Total number of reads: 313280012
Barcodes used: 192

#### "Nov15":
`/illumina_ST-J00101flowcellA_SampleIdNA_RunId0016_LaneId2/Undetermined_S0_L002_R1_001.fastq.gz`
File size: 14,9 GB
Total number of reads: 166348949
Barcodes used: 192

#### "Eu1" ("Eurasia1"):
`/illumina_ST-J00101flowcellA_SampleIdNA_RunId0049_LaneId6/Undetermined_S0_L006_R1_001.fastq.gz`
File size: 17,5 GB
Total number of reads: 200247739 
Barcodes used: 384

#### "Eu2" ("Eurasia2"):
`/illumina_ST-J00101flowcellA_SampleIdNA_RunId0049_LaneId7/Undetermined_S0_L007_R1_001.fastq.gz`
File size: 27,8 GB
Total number of reads: 319991542
Barcodes used: 96


### Reference genome

The *Capsella bursa-pastoris* reference genome used in this project is `Cbp-2-2_contigs.fasta`, provided by Dan Koenig.
A concatenated "pseudogenome" from Cor and Cru referece genomes was also used (`CbpPseudoRef.fasta`). 


### Adaptor barcode sequences

"All-barcodes-all-lanes-8-1-2018.csv" contains a list which adaptors is ligated to which sample in each lane:
- Column 1 is barcode number
- Column 2 is barcode sequence
- following colums refer to each lane, and the rows contain the sample name matching the barcode. If a barcode hasn't been used in lane, the barcode sequence is used again as place holder.

For demultiplexing, each lane has its own txt file (with barcode sequence in line 1 and sample name in line 2):
 - Brcds-Nov14.txt
 - Brcd-Am-Pool1-Switch.txt
 - Brcd-Am-Pool2-Switch.txt
 - Brcds-Nov15.txt
 - Brcd-Eu-Pool1.txt
 - Brcd-Eu-Pool2.txt


## Secondary Analysis

### Step 1: Demultiplexing

Demultiplexing python file (custom): 'Demultiplex-CW-7-4-17.py' (Christina Wesse)
Reads through raw.fastq.gz each read by read, checks for correct formatting (1st line starts with "@", 2nd line is DNA, 3rd line is "+", 4th line is PHRED), looks for the RE cut site, sorts by adaptor barcode and saves reads according to the adaptor sequence, if found in given .txt file.
Generates output as "sample_samplename.fastq.gz" for each barcode including 1 "trash" file ("trash.fastq.gz"; no cut site found) and "no_barcodes.fastq.gz" (cut site but no barcode found). With saving the reads, barcode sequence and corresponding PHRED was cut off, so every read starts with the actual cut site sequence.
Python file also generates *.txt.log with numbers of reads per output fastq.gz file. *: Logfile will be named after input barcodes.txt file for distinct naming.

Execute Demultiplex-CW-7-4-17.py with according barcodes text file for each sequencing lane like this:
- call python script with 2 arguments: 1: path/to/<raw.fastq.gz> 2: path/to/<barcodesfile>.txt
- Barcode input file must contain barcode sequences in column 1 and sample names in column 2 separated by a space
- Example of Barcodes text file:
	GCATT 1583-05-00-00
	TGACT 1198-05-00-00
	GAATG 0966-04-00-00
	TCGAG 1379-10-00-00
	TAGGC 1583-03-00-00
- Default restriction enzyme is KpnI ("REsite", line 73)


The Python script was executed on the HPCC biocluster of the University of California, Riverside, with submitting one Job per lane with `sbatch Job-*name*.sh`:

 - Job-Nov14.sh (see "Nov14"; log file is `Brcds-Nov14.txt.log`)
 
 - Job-Am1-switch.sh (see "America1"; log file is `Brcd-Am-Pool1-Switch.txt.log`)
 
 - Job-Am2-switch.sh (see "America2"; log file is `Brcd-Am-Pool2-Switch.txt.log`)

 - Job-Nov15.sh (see "Nov15"; log file is `Brcds-Nov15.txt.log`)

 - Job-Eu1.sh (see "Eurasia1"; log file is `Brcd-Eu1.txt.log`)
 
 - Job-Eu2.sh (see "Eurasia2"; log file is `Brcd-Eu2.txt.log`)
 
 


## Step 2: Trim with *Trimmomatic* and align to reference genome with *Burrows-Wheeler Aligner* (*BWA*)

To eliminate bad quality reads *Trimmomatic 0.36* (Bolger et al. 2014) was used:

`java -jar trimmomatic-0.36.jar SE ${fileFullPath} -trimlog ${trimLogfile} ${trimmedOutputname} ILLUMINACLIP:${adaptersFile} SLIDINGWINDOW:4:15 MINLEN:75`


The arguments do the following:

 - *SE* for single end data
 - *-trimlog*  generates a logfile with specified name
 - Cut adapter and other illumina-specific sequences with *ILLUMINACLIP* from the read (ILLUMINACLIP:path/to/adapters.fa).
 - Scan the read with a *4*-base wide *sliding window*, cutting when the average quality per base drops below *15*
 - Drop reads below *75* bases long (*MINLEN*)

To align the reads to the reference genome *bwa 0.7.12* was used:
 
BWA (Heng and Durbin 2009) is a software package for mapping sequences against a reference genome. The algorithm BWA-MEM (Max Exact Matches) is designed for Illumina sequence reads from 70bp to 1Mbp. Since *Trimmomatic* was run with `MINLEN:75` (see "Step 2"), the reads should be >=75 bp, so BWA-MEM has been used. It has the best performace for 70-100bp Illumina reads (http://bio-bwa.sourceforge.net/).
Sorting of the alignment file was performed with samtools 1.4.1 (Heng et al. 2009).

`bwa index ${referenceGenome}`
`bwa mem -R "@RG\tID:${file}\tPU:${flowCell}.${lane}.${sampleName}\tPL:ILLUMINA\tSM:${sampleName}" ${referenceGenome} ${trimmedOutputname} > ${sampleName}.sam`
`samtools view -bSu ${sampleName}.sam | \`
`samtools sort - -o ${sampleName}.sorted.bam`
`samtools index ${sampleName}.sorted.bam  ${sampleName}.sorted.bai`

The arguments do the following:

 - *index* the reference genome
 - use the algorithm bwa *mem*
 - *samtools view -bSu*: *b*: Output in the BAM format. *S*: Input is SAM format. u: Output uncompressed BAM. This option saves time spent on compression/decompression and is thus preferred when the output is piped to another samtools command.
 - *sort* the reads and *o*utput to FILE
 - *index* the BAM file


Both, *Trimmomatic* and *BWA* were executed on the biocluster with submitting one Job per lane with `sbatch *name*-xxx-pipeline-trim-bwa.sh`:

Cbp-Reference (Cbp-2-2_contigs.fasta):

 - `Nov14-Job-pipeline-trim-bwa.sh`

 - `Am1-turn-Job-pipeline-trim-bwa.sh` *has not been done with the "switched" samples yet!!

 - `Am2-turn-Job-pipeline-trim-bwa.sh` *has not been done with the "switched" samples yet!!

 - `Nov15-Job-pipeline-trim-bwa.sh`

 - `Eu1-Job-pipeline-trim-bwa.sh` 
 
 - `Eu2-Job-pipeline-trim-bwa.sh`
 
Pseudoreference (CbpPseudoRef.fasta):
 
 - `Job-bwa-PseudoRef-Nov14.sh`
 
 - `Job-bwa-PseudoRef-Am1-switch.sh`
 
 - `Job-bwa-PseudoRef-Am2-switch.sh`
 
 - `Job-bwa-PseudoRef-Nov15.sh`
 
 - `Job-bwa-PseudoRef-Eu1.sh`
 
 - `Job-bwa-PseudoRef-Eu2.sh`
 
 

Log files are `Trimlog_<samplenumber>.txt` for each aligned sample.
Only "real" samples (those demultiplexed fastq.gz with sample names, not barcodes names) were processed. The output is one "samplename".sorted.bam, one "samplename".sorted.bam.bai and one Trimlog_"samplename".txt per sample. 

## Step 3: Call and filter SNPs

The program *freebayes v9.9.2-27-g5d5b8ac* (Garrison and Marth 2012) was used for SNP calling. For better performance the process was parallelized with *freebayes-parallel*:

`bash freebayes-parallel.sh <(python fasta_generate_regions.py ${referenceGenome}.fai 100000) 14 | \`
`--no-indels -f ${referenceGenome} --bam-list ${filesList} | gzip > ${outputname}.gz`

The arguments do the following:

 - *100000* seperates the reference genome into chunks of 100000 bases length to be processed in parallel
 - *14* is the number of parallel processes 
 - *no-indels* ignores indels 
 - *bam-list* reads in a specified txt file (*filesList*) with sample names to process


It was executed with one Job for all data: `sbatch Job-freebayes-alldata-8-6-2017.sh` (deprecated; since there were problems with Cbp-2-2_contigs.fasta)
For pseudoreference: `Job-freebayes-Pseudoref.sh` (with `filesListPseudoRef.txt` as ${filesList})

After that, SNP filtering was performed with *vcftools* version 0.1.13(Danecek et al. 2011). 

Pseudoreference: `Job-filterSNPs-pseudo-8-1-2018.sh`: 
`vcftools --gzvcf ${inputFile} --recode --recode-INFO-all --remove-indels | \`
`--max-missing 0.7 --minQ 30 --min-alleles 2 --max-alleles 2  --out ${outputFile}_1`

The arguments do the following:

 - *recode*: write a new vcf output file
 - *recode-INFO-all*: keep all info from the input file
 - *remove-indels*: remove variants that alter the length of the REF allele
 - *max-missing*: Exclude sites with less than 70% data coverage
 - *minQ*: Includes only sites with Quality value above 30
 - *min-alleles 2 max-alleles 2*: include only biallelic sites



## Tertiary Analysis

to be continued
----------


## References

__Bolger, Anthony M., Marc Lohse, and Bjoern Usadel. "Trimmomatic: a flexible trimmer for Illumina sequence data." Bioinformatics (2014): btu170.

__Li, Heng, and Richard Durbin. "Fast and accurate short read alignment with Burrows–Wheeler transform." Bioinformatics 25.14 (2009): 1754-1760.

__Li, Heng, et al. "The sequence alignment/map format and SAMtools." Bioinformatics 25.16 (2009): 2078-2079.

__Garrison, Erik, and Gabor Marth. "Haplotype-based variant detection from short-read sequencing." arXiv preprint arXiv:1207.3907 (2012).

__Danecek, Petr, et al. "The variant call format and VCFtools." Bioinformatics 27.15 (2011): 2156-2158.
