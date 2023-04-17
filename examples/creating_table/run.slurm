#!/bin/sh
#SBATCH --job-name=TSO_OPA
#SBATCH --nodes=5
#SBATCH --ntasks=160
##SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --partition=debug
#SBATCH --mail-user=eitner@mpia.de  
#SBATCH --mail-type=BEGIN,END,FAIL
    
~/bin/julia create_table.jl