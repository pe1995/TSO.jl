#!/bin/sh
#SBATCH --job-name='TSO'
#SBATCH --nodes=3
#SBATCH --ntasks=100
#SBATCH --time=05:00:00
#SBATCH --mem=92G
    
julia create_table.jl