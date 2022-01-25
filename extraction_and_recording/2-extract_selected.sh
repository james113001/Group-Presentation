#PBS -l walltime=5:00:00
#PBS -l select=1:ncpus=1:mem=20gb
#PBS -N extraction
#PBS -q med-bio

cd /rds/general/user/bbodinie/projects/hda_21-22/live/TDS/General/extraction_and_recoding/scripts
module load anaconda3/personal

ukb_path=/rds/general/user/bbodinie/projects/hda_21-22/live/TDS/General/Data/ukb47946.csv

Rscript 2-extract_selected.R $ukb_path 

