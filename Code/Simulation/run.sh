#!bash

cdir=${PWD}
for init in 0.00  #0.15 0.20 0.25 0.30 #0.05 0.10
do
for eps in 5
do
for rat in 2.9 # 0.5 #0.0 0.5 1.0 2.0 3.0
do
for size in 7.5 # 3.0 3.5 #5.0 #1.5 5
do
for seed in `seq 1 30 `
do
name=init_${init}0_${rat}00_${seed}_${eps}_${size} 
python3 /home/ucapbbm/Scratch/Dumbbells/scripts/generatevesicals/build_LargeVescHorizMass.py 1200000 ${init} ${rat} 1.0 1 ${seed} ${eps} ${size}
python3 /home/ucapbbm/Scratch/Dumbbells/scripts/generatevesicals/gen_submit.py ${cdir} ${name} 
ndir=${cdir}/${name}
name2=Mass30
#rm ${ndir}/submitRe.pbs
#cp -r ${ndir} ${name2}/ 
#rm -r ${ndir}
qsub ${ndir}/submit.pbs
done
done
done
done
done
