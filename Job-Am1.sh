#!/bin/bash -lex
# Exec python script "Demultiplex-CW-7-4-17" with "AmericaPool1" raw fastq.gz and its barcode txt file

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=3-00:00:00     # 3 days
#SBATCH --output=report-Am1-brcds
#SBATCH --mail-user=christina.wesse.osnabrueck@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --job-name="CW-RAD-Am1-brcds"
#SBATCH -p intel # This is the default partition, you can use any of the following; intel, batch, highmem, gpu

newDir="Am1-switch"
rawlane="/rhome/cwesse/shared/OUTSIDE_SEQUENCE_DATA/CHRISTINA_CBP_RAD/illumina_ST-J00101flowcellA_SampleIdNA_RunId0008_LaneId7/Undetermined_S0_L007_R1_001.fastq.gz"
brcdsfile="../Brcd-Am-Pool1-Switch.txt"

mkdir ${newDir}
cd ${newDir}
python ../Demultiplex-CW-7-4-17.py ${rawlane} ${brcdsfile}

echo "done"
