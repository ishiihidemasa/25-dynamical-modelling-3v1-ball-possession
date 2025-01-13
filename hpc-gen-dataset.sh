#!/bin/sh
nohup apptainer exec julia.sif julia --project=. scripts/generate-dataset.jl &> gendata.out