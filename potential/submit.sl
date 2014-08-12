#!/bin/bash

# set the partition where the job will run
#SBATCH --partition=gpu

#SBATCH --ntasks-per-node=4
#SBATCH --ntasks-per-socket=2

# set the number of nodes
#SBATCH --nodes=1

# set the number of GPU cards to use per node
#SBATCH --gres=gpu:1

# set max wallclock time
#SBATCH --time=00:00:10






./runSwitcher 5
