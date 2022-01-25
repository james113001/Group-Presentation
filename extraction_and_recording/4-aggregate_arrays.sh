#PBS -l walltime=5:00:00
#PBS -l select=1:ncpus=1:mem=20gb
#PBS -N aggregating
#PBS -q med-bio

cd /rds/general/user/bbodinie/projects/hda_21-22/live/TDS/General/extraction_and_recoding/scripts
module load anaconda3/personal

Rscript 4-aggregate_arrays.R

# Creating zip file with parameters used for this project
cd ../
project_name=myproject
cp -r parameters parameters_$project_name
zip -r parameters_$project_name.zip parameters_$project_name
rm -rf parameters_$project_name
