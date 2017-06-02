# *Capsella bursa-pastoris* RAD-Seq

Last modified: June 2nd 2017, Christina Wesse

## Primary Analysis
During 6 sessions from Nov 2014 to Dez 2016 tissue from Shepherd's Purse *Capsella bursa-pastoris* (Cbp in the following) was DNA extracted with CTAB and chloroform. Library prep for single end RAD sequencing was performed with 384 individual adaptor barcodes each multiplex. Sequencing was performed in the Genome Center of the Max Planck Institute for Developmental Biology in Tübingen, Germany. DNA extracts and library backup are stored in the Weigel Lab. 

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

#### "Am1-turn"* ("America1"):
`/illumina_ST-J00101flowcellA_SampleIdNA_RunId0008_LaneId7/Undetermined_S0_L007_R1_001.fastq.gz`
File size: 19,1 GB
Total number of reads: 323376947
Barcodes used: 232

#### "Am2-turn"* ("America2"):
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

*: During library prep for America Pools 1 and 2, barcodes plate 1 (1-96) was turned upside down, so the barcodes txt files had to be changed. This is why these data are called "turn". The barcodes txt files (see "All-barcodes-all-lanes.csv") contain the sample names in the corrected order.

### Adaptor barcode sequences

"All-barcodes-all-lanes.csv" contains a list which adaptors is ligated to which sample in each lane:
- Column 1 is barcode number
- Column 2 is barcode sequence
- following colums refer to each lane, and the rows contain the sample name matching the barcode. If a barcode hasn't been used in lane, the barcode sequence is used again as place holder.

For demultiplexing, each lane has its own txt file (with barcode sequence in line 1 and sample name in line 2):
 - Brcds-Nov14.txt
 - Brcd-Am-Pool1-turned-plate-1.txt
 - Brcd-Am-Pool2-turned-plate-1.txt
 - Brcds-Nov15.txt
 - Brcd-Eu-Pool1.txt
 - Brcd-Eu-Pool2.txt

## Secondary Analysis

### Step 1: Demultiplexing

Demultiplexing python file (custom): 'Demultiplex-CW-7-4-17.py' (Christina Wesse)
Reads through raw.fastq.gz each read by read, checks for correct formatting (1st line starts with "@", 2nd line is DNA, 3rd line is "+", 4th line is PHRED), looks for the RE cut site, sorts by adaptor barcode and saves reads according to the adaptor sequence, if found in given .txt file.
Generates output as "sample_samplename.fastq.gz" for each barcode including 1 "trash" file ("trash.fastq.gz"; no cut site found) and "no_barcodes.fastq.gz" (cut site but no barcode found). With saving the reads, barcode sequence and corresponding PHRED was cut off, so every read starts with the actual cut site sequence.
Python file also generates *.txt.log with numbers of reads per output fastq.gz file. *: Logfile will be named after input barcodes.txt file for distinct naming.

The Python script was executed on the biocluster with submitting one Job per lane with `sbatch Job-*name*.sh`:

 - Job-Nov14.sh (see "Nov14"; log file is `Brcds-Nov14.txt.log`)

 - Job-Am1.sh (see "America1"; log file is `Brcd-Am1.txt.log`)

 - Job-Am2.sh (see "America2"; log file is `Brcd-Am1.txt.log`)

 - Job-Nov15.sh (see "Nov15"; log file is `Brcds-Nov15.txt.log`)

 - Job-Eu1.sh (see "Eurasia1"; log file is `Brcd-Eu1.txt.log`)
 
 - Job-Eu2.sh (see "Eurasia2"; log file is `Brcd-Eu2.txt.log`)


## Step 2: Trim with Trimmomatic and align to reference genome with Burrows-Wheeler Aligner (BWA)

To eliminate bad quality reads *Trimmomatic 0.36* was used:

`java -jar trimmomatic-0.36.jar SE ${fileFullPath} -trimlog ${trimLogfile} ${trimmedOutputname} ILLUMINACLIP:${adaptersFile} SLIDINGWINDOW:4:15 MINLEN:75`


The arguments do the following:

 - SE for single end data
 - -trimlog  generated a logfile with specified name
 - Cut adapter and other illumina-specific sequences from the read (ILLUMINACLIP:path/to/adapters.fa).
 - Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15 (SLIDINGWINDOW:4:15)
 - Drop reads below the 75 bases long (MINLEN:75)

To align the reads to the reference genome *bwa 0.7.12* was used:

BWA is a software package for mapping sequences against a reference genome. The algorithm BWA-MEM (Max Exact Matches) is designed for Illumina sequence reads from 70bp to 1Mbp. Since *Trimmomatic* was run with `MINLEN:75` (see "Step 2"), the reads should be >=75 bp, so BWA-MEM has been used. It has the best performace for 70-100bp Illumina reads (http://bio-bwa.sourceforge.net/).

`bwa index ${referenceGenome}`
`bwa mem -R "@RG\tID:${file}\tPU:${flowCell}.${lane}.${sampleName}\tPL:ILLUMINA\tSM:${sampleName}" ${referenceGenome} ${trimmedOutputname} > ${sampleName}.sam`
`samtools view -bSu ${sampleName}.sam | \`
`samtools sort - -o ${sampleName}.sorted.bam`
`samtools index ${sampleName}.sorted.bam  ${sampleName}.sorted.bai`

The programs were executed on the biocluster with submitting one Job per lane with `sbatch *name*-Job-pipeline-trim-bwa.sh`:

 - `Nov14-Job-pipeline-trim-bwa.sh`

 - `Am1-turn-Job-pipeline-trim-bwa.sh`

 - `Am2-turn-Job-pipeline-trim-bwa.sh` 

 - `Nov15-Job-pipeline-trim-bwa.sh`

 - `Eu1-Job-pipeline-trim-bwa.sh` 
 
 - `Eu2-Job-pipeline-trim-bwa.sh`

Only "real" samples (those demultiplexed fastq.gz with sample names, not barcodes names) were processed. The output is one "samplename".sorted.bam, one "samplename".sorted.bam.bai and one Trimlog_"samplename".txt per sample. 

## Step 3: Call and filter SNPs

IN PROGRESS

## Tertiary Analysis

----------


## References

__Bolger, Anthony M., Marc Lohse, and Bjoern Usadel. "Trimmomatic: a flexible trimmer for Illumina sequence data." Bioinformatics (2014): btu170.

__Li, Heng, and Richard Durbin. "Fast and accurate short read alignment with Burrows–Wheeler transform." Bioinformatics 25.14 (2009): 1754-1760.

__Li, Heng, et al. "The sequence alignment/map format and SAMtools." Bioinformatics 25.16 (2009): 2078-2079.

__Garrison, Erik, and Gabor Marth. "Haplotype-based variant detection from short-read sequencing." arXiv preprint arXiv:1207.3907 (2012).

__Danecek, Petr, et al. "The variant call format and VCFtools." Bioinformatics 27.15 (2011): 2156-2158.
