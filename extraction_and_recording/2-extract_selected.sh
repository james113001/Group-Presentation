#PBS -l walltime=5:00:00
#PBS -l select=1:ncpus=1:mem=40gb
#PBS -N extraction
#PBS -q med-bio

cd /rds/general/project/hda_21-22/live/TDS/Group_8/Group-Presentation/extraction_and_recording/scripts
module load anaconda3/personal

ukb_path=/rds/general/project/hda_21-22/live/TDS/General/Data/ukb47946.csv

Rscript 2-extract_selected.R $ukb_path 

