#!/bin/bash -l
# Exec python script "Demultiplex-CW-7-4-17" with "AmericaPool1" raw fastq.gz and its barcode txt file

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=3-00:00:00     # 3 days
#SBATCH --output=report-Am1
#SBATCH --mail-user=christina.wesse.osnabrueck@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --job-name="CW-RAD-Am1"
#SBATCH -p intel # This is the default partition, you can use any of the following; intel, batch, highmem, gpu

cd /bigdata/koeniglab/cwesse/CBP-RAD-CW/Am1-turn
python ../Demultiplex-CW-7-4-17.py /rhome/cwesse/shared/OUTSIDE_SEQUENCE_DATA/CHRISTINA_CBP_RAD/illumina_ST-J00101flowcellA_SampleIdNA_RunId0008_LaneId7/Undetermined_S0_L007_R1_001.fastq.gz /bigdata/koeniglab/cwesse/CBP-RAD-CW/Brcd-Am-Pool1-turned-plate-1.txt

