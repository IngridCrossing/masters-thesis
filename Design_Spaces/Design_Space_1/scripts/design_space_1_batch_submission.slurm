#!/bin/bash

for survey_designs in {1..10}
do

  for det_prob_options in {1..10}
  do
    
      sbatch Design_Spaces/Design_Space_1/scripts/design_space_1_job_submission.slurm $survey_designs $det_prob_options
      
    done
done