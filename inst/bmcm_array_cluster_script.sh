#!/bin/bash -l

# Example jobscript to run a serial array job

# Request ten minutes of wallclock time (format hours:minutes:seconds)
#$ -l h_rt=2:10:0

# Request 1 gigabyte of RAM (must be an integer followed by M, G, or T)
#$ -l mem=40G

# Request 15 gigabyte of TMPDIR space (default is 10 GB - remove if cluster is diskless)
#$ -l tmpfs=15G

# Set up the job array. In this instance we have requested 1000 tasks
# numbered 1 to 1000
#$ -t 1-2

# Set the name of the job
#$ -N R_job_bmcm

# Set the working directory to somewhere in your scratch space
#  This is a necessary step as compute nodes cannot write to $HOME
# Replace "<your_UCL_id>" with your UCL user ID
#$ -wd /home/sejjng1/Scratch

## Your work should be done in $TMPDIR 
#cd $TMPDIR

# Load the R module and run the R program
module -f unload compilers mpi gcc-libs
module load r/4.4.0-openblas/gnu-10.2.0

# add library folder to path
export R_LIBS=/lustre/home/sejjng1/R/x86_64-pc-linux-gnu-library/4.2:$R_LIBS  

# Get the task ID
TASK_ID=$SGE_TASK_ID

# select for hierarchical or separate models:
#R --no-save --args ${TASK_ID} ${SCENARIO_ID} < /home/sejjng1/Scratch/bmcm/cluster_script.R
R --no-save --args ${TASK_ID} ${SCENARIO_ID} < /home/sejjng1/Scratch/bmcm/cluster_script_separate.R

# Preferably, tar-up (archive) all output files onto the shared scratch area
#tar -zcvf $HOME/Scratch/R_output/files_from_job_$JOB_ID.tar.gz $TMPDIR
