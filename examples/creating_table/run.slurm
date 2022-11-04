#!/bin/sh
#SBATCH --job-name='TSO'
#SBATCH --nodes=4
#SBATCH --ntasks=128
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --partition=debug
    
julia create_table.jl