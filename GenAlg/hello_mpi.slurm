#! /bin/bash
### Start parallel python jobs

### Cluster-Job Configuration
#SBATCH --cluster=omics
#SBATCH --partition=shortterm
#SBATCH --mail-type=ALL
#SBATCH --mail-user=michael.olbrich@student.uni-luebeck.de
#SBATCH --job-name=hello_mpi

#SBATCH --time=0:02:00
#SBATCH --mem-per-cpu=1G
#SBATCH --ntasks=4

#SBATCH --no-requeue

source activate ga_workspace

mpirun python hello_mpi.py