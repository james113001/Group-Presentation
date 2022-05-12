#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=1:mem=10gb
#PBS -N sevnodes1core
#PBS -J 1-8


cd /rds/general/project/hda_21-22/live/TDS/Group_8/Group-Presentation/extraction_and_recording/'ANALYSIS - descriptive'

module load anaconda3/personal

ichunk=$PBS_ARRAY_INDEX
nchunks=8

Rscript 'Summary file [1:30]'.R $nchunks  $ichunk