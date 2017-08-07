#!/bin/bash
#SBATCH --job-name=rot_MPtest
#SBATCH --account=commons
#SBATCH --partition=commons
#SBATCH --nodes=1
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=500M
#SBATCH --time=00:30:00
#SBATCH --mail-user=sk71@rice.edu
#SBATCH --mail-type=ALL
 
echo "My job ran on:"
echo $SLURM_NODELIST
    cd $HOME/BDMP
    srun ./b2 ./input/sample_srfrpl_openMP ./output/output_BDMP_test.dat