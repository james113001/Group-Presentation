#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=1:mem=10gb
#PBS -N plot
#PBS -q med-bio


cd /rds/general/project/hda_21-22/live/TDS/Group_8/Group-Presentation/extraction_and_recording/'ANALYSIS - descriptive'

module load anaconda3/personal

Rscript plots_pvals.R