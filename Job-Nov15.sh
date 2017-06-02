#!/bin/bash -l
# Exec python script "Demultiplex-CW-7-4-17" with "Nov15" raw fastq.gz and its barcode txt file

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=3-00:00:00     # 3 days
#SBATCH --output=report-Nov15
#SBATCH --mail-user=christina.wesse.osnabrueck@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --job-name="CW-RAD-Nov15"
#SBATCH -p intel # This is the default partition, you can use any of the following; intel, batch, highmem, gpu

cd /bigdata/koeniglab/cwesse/CBP-RAD-CW/Nov15
python ../Demultiplex-CW-7-4-17.py /rhome/cwesse/shared/OUTSIDE_SEQUENCE_DATA/CHRISTINA_CBP_RAD/illumina_ST-J00101flowcellA_SampleIdNA_RunId0016_LaneId2/Undetermined_S0_L002_R1_001.fastq.gz /bigdata/koeniglab/cwesse/CBP-RAD-CW/Brcds-Nov15.txt

