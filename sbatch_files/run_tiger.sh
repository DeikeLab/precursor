#!/bin/bash
#SBATCH --job-name=720_0p3          # create a short name for your job
#SBATCH --nodes=4                   # node count
#SBATCH --ntasks-per-node=40        # number of tasks per node
#SBATCH --cpus-per-task=1           # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G            # memory per cpu-core (4G is default)
#SBATCH --time=23:59:00             # total run time limit (HH:MM:SS)
##SBATCH --time=00:30:00             # total run time limit (HH:MM:SS)
####SBATCH --mail-type=begin        # send email when job begins
####SBATCH --mail-type=end          # send email when job ends
####SBATCH --mail-user=ns8802@princeton.edu
#
module load openmpi/gcc/3.1.5/64
#
ulimit -s unlimited
#
# Input parameters
#
Re_ast=720.0;                 
ak=0.3;                       
r_L0lam=4.0;                       
rho_r=0.001225;           
mu_r=2.2471881948940954e-06; 
MAXLEVEL=9;                   
MINLEVEL=6;                  
f_uemax=0.3;                 
Tf_=3.14;                    
do_en=0;                      
st_wave=1;                    
rand_num=2;                     
dump_now=0;                   
N_mod=192;                     
#
srun precursor -m 23:59:00 $Re_ast $ak $r_L0lam $rho_r $mu_r $MAXLEVEL $MINLEVEL $f_uemax \
	                   $Tf_ $do_en $st_wave $rand_num $dump_now $N_mod > out.log 2>&1
#
