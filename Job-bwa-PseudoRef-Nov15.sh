#!/bin/bash -lex
# Exec Trimmomatic and align trimmed reads with bwa
# the CbpPseudoRef.fasta was concatenated from Co_contigs.fasta and Crubella_183_v1.fa

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=3G
#SBATCH --time=5-00:00:00     # 5 days
#SBATCH --output=report-bwa-pseudoref-Nov15
#SBATCH --mail-user=christina.wesse.osnabrueck@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --job-name="pseudoref-Nov15-turn-CBP-Trim-BWA"
#SBATCH -p intel # This is the default partition, you can use any of the following; intel, batch, highmem, gpu

module load bwa
module load samtools

adaptersFile="TruSeq3-SE.fa"
referenceGenome="CbpPseudoRef.fasta"
samFolder="pseudoref-samresults-Nov15"

cd ..
mkdir -p ${samFolder}
cp -u refGenomes/${referenceGenome} ${samFolder}

# Index the reference genome

cd ${samFolder}
if  [ ! -e "CbpPseudoRef.fasta" ] ||  
	[ ! -e "CbpPseudoRef.fasta.ann" ] ||  
	[ ! -e "CbpPseudoRef.fasta.pac" ] ||  
	[ ! -e "CCbpPseudoRef.fasta.amb" ] ||  
	[ ! -e "CbpPseudoRef.fasta.bwt" ] ||  
	[ ! -e "CbpPseudoRef.fasta.sa" ] ||
	[ ! -e "CbpPseudoRef.fasta.fai" ]; then

	echo "Generating files from ${referenceGenome}"
	bwa index ${referenceGenome}
	samtools faidx ${referenceGenome}
	echo "Generating files from ${referenceGenome} done"
else
	echo "Skipping Generating files from ${referenceGenome}"
fi
cd ..

# $1 first argument is folder
# Trimmomatic bwa pipeline
create_sam_file()
{
  files=$(find $1 | grep "sample" | grep -v "T" | grep -v "G" | grep -v "C")
  numberOfSamples=$(ls ${files} | wc -l)
  echo "Number of samples in folder ${1} is ${numberOfSamples}"
  pushd ${samFolder}
	for file in ${files}; do
		
		fileFullPath="../${file}"										   # example: ../demux-jobs/Nov15/sample_1886-02-00-00.fastq.gz
		fileName=$(echo ${fileFullPath} | cut -d '/' -f4  ) 			   # example: sample_1886-02-00-00.fastq.gz
		sampleName=$(echo ${fileName} | sed 's#sample_##g' | sed 's#.fastq.gz##g' )  # example: 1886-02-00-00
		trimLogfile="Trimlog_${sampleName}.txt"							   # example: Trimlog_1886-02-00-00.txt
		readHeader=$(zcat ${fileFullPath} | grep -m1 "@HISEQ\|@ST")		   # get the read header from one read from the sample file; example: @ST-J00101:16:H35T2BBXX:2:1101:11606:1279 1:N:0:0
		flowCell=$(echo ${readHeader} | cut -d ':' -f2)					   # take the colon as separator and take the second element from string, which corresponds to the flow cell ID. example: 16
		lane=$(echo ${readHeader} | cut -d ':' -f4)						   # as above, but take the 4th element, which corresponds to the lane of the flow cell. example: 2
		trimmedOutputname=$(echo ${fileName} | sed 's#sample#filtered#g')  # example: filtered_1886-02-00-00.fastq.gz
		
		echo "Taking file ${fileFullPath}"
		echo "We are in $(pwd)"
		
		echo "execute Trimmomatic"
		java -jar ../align-jobs/trimmomatic-0.36.jar \
		  SE ${fileFullPath} \
		  -trimlog ${trimLogfile} \
		  ${trimmedOutputname} \
		  ILLUMINACLIP:"../align-jobs/${adaptersFile}":2:30:10 \
		  SLIDINGWINDOW:4:15 \
		  MINLEN:75
		echo "Trimmomatic done"
		echo "Start alignin with bwa" 
		bwa mem -R "@RG\tID:${file}\tPU:${flowCell}.${lane}.${sampleName}\tPL:ILLUMINA\tSM:${sampleName}" ${referenceGenome} ${trimmedOutputname} > ${sampleName}.sam  # creates sam file
		echo "successfully created sam file. Start creating sorted.bam"
		samtools view -bSu ${sampleName}.sam | \
		samtools sort - -o ${sampleName}.sorted.bam 				 	# takes the piped bam file and sorts it
		rm ${trimmedOutputname}		
		rm ${sampleName}.sam
		echo "Creating sorted bam file done"
				
		echo "index the sorted bam file"		
		samtools index ${sampleName}.sorted.bam  ${sampleName}.sorted.bai
		echo "Creating bam index done"

		echo "Processing file ${fileFullPath} done"
	done 
  popd
  echo "create data of file $1 done" 
}

#create_sam_file "demux-jobs/Am1-turn/"
#create_sam_file "demux-jobs/Am2-turn/"
#create_sam_file "demux-jobs/Eu1/"
#create_sam_file "demux-jobs/Eu2/"
#create_sam_file "demux-jobs/Nov14/"
create_sam_file "demux-jobs/Nov15/"
