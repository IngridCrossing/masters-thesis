#!/bin/bash
# Created by the University of Melbourne job script generator for SLURM
# Tue Jul 12 2022 18:58:29 GMT+1000 (Australian Eastern Standard Time)

# Partition for the job:
#SBATCH --partition=physical

# Multithreaded (SMP) job: must run on one node 
#SBATCH --nodes=1

# The name of the job:
#SBATCH --job-name="design_space_3"

# The project ID which this job should run under:
#SBATCH --account="punim1753" #Choose project

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3

# Use this email address:
#SBATCH --mail-user= #Insert email here

# Send yourself an email when the job:
# aborts abnormally (fails)
#SBATCH --mail-type=FAIL
# begins
#SBATCH --mail-type=BEGIN
# ends successfully
#SBATCH --mail-type=END

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=3-0:0:00

# check that the script is launched with sbatch
if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch"
   exit 1
fi

# Run the job from this directory:
#cd #uncomment and Insert directory here


survey_designs=$1
det_prob_options=$2


# The modules to load:
module load foss/2020b
module load r/4.1.0
module load jags/4.3.0

# The job command(s):
Rscript Design_Spaces/Design_Space_3/scripts/run_design_space_3.R $survey_designs $det_prob_options

##DO NOT ADD/EDIT BEYOND THIS LINE##
##Job monitor command to list the resource usage
my-job-stats -a -n -s
