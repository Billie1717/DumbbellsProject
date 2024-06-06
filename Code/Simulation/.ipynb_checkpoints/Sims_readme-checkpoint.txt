#------------------------------------------------------------------------------#
Simulations "protocol"
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
Example input and run : /run_example/
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
Python script to generate inputs (horizontal dumbbell) to create particles.in and in.local:

python3 /Users/billiemeadowcroft/Dropbox/PhDGithub/RemoteDumbbells/Code/generatevesicles/build_LargeVescHorizMass.py 1200000 ${init} ${rat} 1.0 1 ${seed} ${eps} ${size}
eg:
python3 /Users/billiemeadowcroft/Dropbox/PhDGithub/RemoteDumbbells/Code/build_LargeVescHorizMass.py 1200000 0.00 2.9 1.0 1 1 5.0 7.5
(same but for build_LargeVescAndMass.py for vericle dumbbell starting position)
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
LAMMPS version I used: 
http://lammps.sandia.gov/tars/lammps-16Feb16.tar.gz
gerun /home/ucapbbm/Scratch/np-toolkit/lammps/src/lmp_mpi -in %s/in.local
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
Parameter sweeps 


1. Many seeds with no solute, and the above parameters (Dumbbell-Membrane epsilon = 5, Dumbbell Size = 7.5, Dumbbell Mass = 30) to get pathways data
-- THESE RUNS ARE STORED IN FOLDERS *** .. in myriad?
A) /home/ucapbbm/Scratch/Dumbbells/scripts/runs/LargeVesc2/VerticleAttempt/Mass30
B) /home/ucapbbm/Scratch/Dumbbells/scripts/runs/LargeVesc2/HorizontalWtMass

2. Solute varied (for tension data) with init. = 0.050 
rat : 0.250 0.500 0.750 1.000 1.250 1.750 2.000 2.250 2.750 3.000
Did this for 2 sizes (& 2 epsilons) :
Size = {5.0 7.5} --> Eps = {8 5}

-- THESE RUNS ARE STORED IN FOLDER ***
A) /home/ucapbbm/Scratch/Dumbbells/scripts/runs/Mass30
