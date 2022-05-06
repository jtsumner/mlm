#! /bin/bash

slurm_out="slurm_mlm.out"

if [ $slurm_out ]
then
    cat $slurm_out >> results/slurm_mlm_old.out
fi