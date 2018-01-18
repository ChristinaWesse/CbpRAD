#!/bin/bash -lex
# Variant calling with freebayes

#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=20G
#SBATCH --time=7-00:00:00     # 7 days
#SBATCH --output=report-pseudoref-snpcall
#SBATCH --mail-user=christina.wesse.osnabrueck@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --job-name="pseudoref-snpcall"
#SBATCH -p intel # This is the default partition, you can use any of the following; intel, batch, highmem, gpu

module load freebayes
module load parallel
module load vcflib

referenceGenome="../refGenomes/CbpPseudoRef.fasta"
echo ${referenceGenome}

filesList="filesListPseudoRef.txt"
cat ${filesList}

outDate=$(date +"%d-%m-%y")
outTime=$(date +"%T")
outputname="pseudoref-${outDate}-${outTime}.vcf"

bash freebayes-parallel.sh <(python fasta_generate_regions.py ${referenceGenome}.fai 100000) 14 --no-indels -f ${referenceGenome} --bam-list ${filesList} | gzip > ${outputname}.gz

echo "process done"

