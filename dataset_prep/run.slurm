#!/bin/bash
#SBATCH --job-name=zip
#SBATCH --partition=56c1t8n
#SBATCH --ntasks-per-node=16
#SBATCH --account=zhangqunjie_lab
#SBATCH --error=%j.error.log
#SBATCH --output=%j.stdout.log
source ~/.bashrc
echo "----------------------------------------------"
echo "|                 [SUMMARY]                   |"
echo "----------------------------------------------"
echo "Job Name: " $SLURM_JOB_NAME
echo "Job ID: " $SLURM_JOBID
echo "User: " $USER
echo "Account: " $SLURM_JOB_ACCOUNT
echo "Partition: " $SLURM_JOB_PARTITION
echo "Node: " $SLURM_NODELIST
echo "Qos: " $SLURM_JOB_QOS
echo "PWD: " $PWD 
echo "----------------------------------------------"
echo "[INFO] Start - "`date`
echo "----------------------------------------------"

# module load compiler/gcc/13.3.0
# conda activate test
# Run ..
python HaplotypeUtils.py

echo "----------------------------------------------"
echo "[INFO] End - "`date`
echo "Done."

# Ver.20240902 BLOSCAU