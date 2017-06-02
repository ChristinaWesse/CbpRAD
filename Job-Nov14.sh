#!/bin/bash -l
# Exec python script "Demultiplex-CW-7-4-17" with "Nov14" raw fastq.gz and its barcode txt file

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=2-00:00:00     # 2 days
#SBATCH --output=report-Nov14
#SBATCH --mail-user=christina.wesse.osnabrueck@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --job-name="CW-RAD-Nov14"
#SBATCH -p intel # This is the default partition, you can use any of the following; intel, batch, highmem, gpu

cd /bigdata/koeniglab/cwesse/CBP-RAD-CW/Nov14
python ../Demultiplex-CW-7-4-17.py /rhome/cwesse/shared/OUTSIDE_SEQUENCE_DATA/CHRISTINA_CBP_RAD/illumina_SN7001143flowcellB_SampleIdNA_RunId0213_LaneId5/lane5_NoIndex_L005_R1_001.fastq.gz /bigdata/koeniglab/cwesse/CBP-RAD-CW/Brcds-Nov14.txt

