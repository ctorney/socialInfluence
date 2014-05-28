
#!/bin/sh

current_time=$(date "+%Y.%m.%d-%H.%M")
for i in `seq 0 5`
do
cat <<EOS | qsub -
#!/bin/sh

#PBS -l nodes=1:ppn=4:gpus=1,walltime=100:00:00
#PBS -j oe
#PBS -q gpu
module load cudatoolkit/5.5.22
module load intel-mkl/11.0/1/64
export CUDA_VISIBLE_DEVICES=\`grep \$HOSTNAME \$PBS_GPUFILE | awk -Fu '{printf A\$2;A=","}'\`

cd \$PBS_O_WORKDIR
./runSwitcher $i
EOS
done
