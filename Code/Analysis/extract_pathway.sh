#!/bin/bash

cdir=${PWD}
for rat in 2.9 
do
for init in 0.00
do
for Eps in 5
do
for size in 7.5
do
echo ${Eps} ${init}
for seed in `seq 1 30`
do
name=init_${init}0_${rat}00_${seed}_${Eps}_${size} 
#for eps in 6
#do
#for R in 30
#do 
echo ${seed}
ndir=${cdir}/${name}
cd ${ndir}
python3 /home/ucapbbm/Scratch/Dumbbells/scripts/analysis/dumbbellpathway.py output.xyz 1.0 
cp ${ndir}/theta_d.dat ${cdir}/pathwaysMass30LargeVesc/theta_d_${init}_${rat}_${seed}_${Eps}_${size}.dat
#rm ${cdir}/pathwaysMass30LargeVesc/theta_d_${init}_${rat}_${seed}_${Eps}_${size}_re.dat
#python3 /home/ucapbbm/Scratch/Dumbbells/scripts/analysis/dumbbellpathway.py output2.xyz 1.0
#cp ${ndir}/theta_d.dat ${cdir}/pathwaysHorizMass30LargeVesc/theta_d_${init}_${rat}_${seed}_${Eps}_${size}_re2.dat
done
done
done
done
done
#done
