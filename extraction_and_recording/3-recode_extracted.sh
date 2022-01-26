#PBS -l walltime=5:00:00
#PBS -l select=1:ncpus=1:mem=20gb
#PBS -N recoding
#PBS -q med-bio

cd /rds/general/project/hda_21-22/live/TDS/Group_8/Group-Presentation/extraction_and_recording/scripts
module load anaconda3/personal

Rscript 3-recode_variables.R

