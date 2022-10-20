#!/bin/bash

for survey_designs in {1..20}
do

  for det_prob_options in {1..10}
  do
    
      sbatch Design/Spaces/Design_Space_2/scripts/design_space_2_job_submission.slurm $survey_designs $det_prob_options
      
    done
done