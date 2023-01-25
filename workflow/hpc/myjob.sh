#!/bin/bash 
# properties = {properties}

source /home/univ4859/.bashrc

source /etc/profile.d/modules.sh
module purge
module load apps/amber/18/GPU-CUDA10

{exec_job}