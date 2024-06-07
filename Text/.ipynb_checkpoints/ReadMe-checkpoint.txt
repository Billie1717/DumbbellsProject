After a lot of exploring of different params:
size of vesical, mass of dumbbell: mass of lipid ratio, Dumbbell-membrane interaction strength (eps), Dumbbell size, ratio of solute in:out

We found a set of params whose pathway data matches that of experiment

IE: 
Size of dumbbell : 7.5 (=sigma)
Mass of dumbbell : 30
Epsilon : 5 (i think)
Size of vesical : 9722 particles (not sure what size this is... ) R= 25.331
Solute : None

In myriad, the horizontally starting dumbbell data seems to be in folder:
/home/ucapbbm/Scratch/Dumbbells/scripts/runs/LargeVesc2/HorizontalWtMass

I think the vericle data is in:
/home/ucapbbm/Scratch/Dumbbells/scripts/runs/LargeVesc2/VerticleAttempt/Mass30
eg :
sftp://ucapbbm@myriad.rc.ucl.ac.uk/home/ucapbbm/Scratch/Dumbbells/scripts/runs/LargeVesc2/VerticleAttempt/Mass30/init_0.000_3.000_1_5_7.5/output1.xyz

Analysis code is in:
/home/ucapbbm/Scratch/Dumbbells/scripts/analysis/extract_pathway.sh


But I should check the params in /home/ucapbbm/Scratch/Dumbbells/scripts/runs/Mass30 as this was the last directory to be mentioned in Build_LargeVescHorizwtMass.py


TO DO:
    -Analysis of verticle runs (check its run long enough)
    -make new plotting script for specific 'pathways' EG 1) Seed 14 horiz->lobe 2) Seed 1 horiz -> sink (with time cut off) (save these .xyz runs on computer too)
    -Start creating paper figure
    

access on Ovito to Myriad:

sftp://ucapbbm@myriad.rc.ucl.ac.uk/home/ucapbbm/Scratch/Dumbbells/scripts/runs/LargeVesc2/HorizontalWtMass/init_0.000_2.900_26_5_7.5/output.xyz
UPDATED:
sftp://ucapbbm@myriad.rc.ucl.ac.uk/home/ucapbbm/Scratch/Dumbbells/scripts/runs/UsedForSnapshots/Horizontal/init_0.000_2.900_26_5_7.5/output.xyz

Screenshots from:

#-------- Horizontal -----------#

Lobe End: sftp://ucapbbm@myriad.rc.ucl.ac.uk/home/ucapbbm/Scratch/Dumbbells/scripts/runs/LargeVesc2/HorizontalWtMass/init_0.000_3.000_3_5_7.5/output.xyz
UPDATED :
sftp://ucapbbm@myriad.rc.ucl.ac.uk/home/ucapbbm/Scratch/Dumbbells/scripts/runs/UsedForSnapshots/Horizontal/init_0.000_3.000_3_5_7.5/output.xyz

Sink End: sftp://ucapbbm@myriad.rc.ucl.ac.uk/home/ucapbbm/Scratch/Dumbbells/scripts/runs/LargeVesc2/HorizontalWtMass/init_0.000_3.000_10_5_7.5/output.xyz
UPDATED:
sftp://ucapbbm@myriad.rc.ucl.ac.uk/home/ucapbbm/Scratch/Dumbbells/scripts/runs/UsedForSnapshots/Horizontal/init_0.000_3.000_10_5_7.5/output.xyz


#-------- Verticle -----------#
/home/ucapbbm/Scratch/Dumbbells/scripts/runs/LargeVesc2/VerticleAttempt/Mass30/init_0.000_3.000_3_5_7.5
sftp://ucapbbm@myriad.rc.ucl.ac.uk/home/ucapbbm/Scratch/Dumbbells/scripts/runs/UsedForSnapshots/Vericle/init_0.000_3.000_3_5_7.5/output1.xyz


#---------- Solute -----------#

EG :InitRat/init_0.150_0.500_21_9_5.0
EG :InitRat/init_0.150_3.000_41_9_5.0
init_${init}0_${rat}00_${seed}_${eps}_${size}

NOTE :
All .xyz files have been removed in /home/ucapbbm/Scratch/Dumbbells/scripts/runs/InitRat/ apart from seeds 6: init_*_*_6_*_*
All .xyz files have been removed in /home/ucapbbm/Scratch/Dumbbells/scripts/runs/Mass30/ apart from seeds 6: init_*_*_6_*_*
All .xyz files have been removed in /home/ucapbbm/Scratch/Dumbbells/scripts/runs/LargeVesc/ apart from seeds 6: init_*_*_6_*_*
All .xyz files have been removed in /home/ucapbbm/Scratch/Dumbbells/scripts/runs/LargeVesc2/ apart from seeds 6: init_*_*_6_*_*
(NOT ANY OTHER FOLDER IN /LargeVesc2 or any other folder in /runs --> Do this when you have collected dumbbells data)


MY USEFUL PAPER DATA:

Seed 14, 1 lobe from horizontal starting:  
Documents/PHD/Dumbbells/Data/pathwaysHorizMass30LargeVesc/theta_d_0.00_2.9_14_5_7.5_re.dat
Documents/PHD/Dumbbells/Data/pathwaysHorizMass30LargeVesc/theta_d_0.00_2.9_14_5_7.5_re2.dat

From myriad data: /home/ucapbbm/Scratch/Dumbbells/scripts/runs/LargeVesc2/HorizontalWtMass/inits with rat 2.9)

state E = HorizontalWtMass/init_0.000_3.000_3_5_7.5/output2.xyz

states A,B,G = VerticleAttempt/Mass30/init_0.000_3.000_3_5_7.5/output1.xyz

states C D H J K = init_0.000_2.900_26_5_7.5/output1.xyz & output2.xyz







(Cisco anyconnect -> myriad -> (DONT HAVE A PROJECT CODE. -> NOTHING IS STORED HERE)
To access simulations, see instructions for logging into RDSS (from Myriad)
https://www.ucl.ac.uk/isd/how-to/rdss-myriad-data-storage-transfer-service
More info about storage:
https://www.ucl.ac.uk/isd/live-storage-access-guide
Info about data storage:
https://www.rc.ucl.ac.uk/docs/