#!/bin/bash
cdir=${PWD}
#if [ -d  wrappingR30N ]; then
#    rm -Rf wrappingR30N
#fi
#mkdir wrappingR20_25N
 

for init in 0.050 #0.000 0.050 0.100 0.150 0.200 0.250 0.300 
do
for seed in `seq 1 30`
do
for rat in 0.250 0.500 0.750 1.000 1.250 1.500 1.750 2.000 2.250 2.500 2.750 3.000 #0.000 0.250 0.500 0.750 1.000 1.250 1.750 2.250 2.750 
do
for R in 5.0 #5.0 5.5 6.0 6.5 #7.0 7.5
do
for eps in 8 #6 7 8
do
echo ${R} ${seed}
ndir=init_${init}_${rat}_${seed}_${eps}_${R}/
cd ${ndir} 
#cmd=python3 /home/ucapbbm/Scratch/Dumbbells/scripts/analysis/wrapping.py output.xyz wrapping_${eps}_${R}_init_${init}_rat_${rat}_seed_${seed}.dat
#${cmd} ${bc}
python3 /home/ucapbbm/Scratch/Dumbbells/scripts/analysis/wrapping.py output2.xyz wrapping_${eps}_${R}_init_${init}_rat_${rat}_seed_${seed}.dat ${R}
#python3 /home/ucapbbm/Scratch/Dumbbells/scripts/analysis/wrapping.py output.xyz wrapping_${eps}_${R}_init_${init}_rat_${rat}_seed_${seed}.dat  ${bc}
cp wrapping_${eps}_${R}_init_${init}_rat_${rat}_seed_${seed}.dat ${cdir}/WrappingInitRatMass30_2
#mv ${cdir}/WrappingEpsSizeMass30/wrapping_${eps}_${R}_init_${init}_rat_${rat}_seed_${seed}.dat ${cdir}/WrappingInitRatMass30
cd ${cdir}
done 
done
done
done
done
