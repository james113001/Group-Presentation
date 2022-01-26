#PBS -l walltime=5:00:00
#PBS -l select=1:ncpus=1:mem=20gb
#PBS -N dict
pwd
cd /rds/general/project/hda_21-22/live/TDS/Group_8/Group-Presentation/extraction_and_recording/scripts
module load anaconda3/personal

ukb_path=/rds/general/project/hda_21-22/live/TDS/General/Data/ukb47946.csv

Rscript 1-make_data_dict.R $ukb_path

