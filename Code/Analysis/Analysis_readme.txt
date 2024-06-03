#------------------------------------------------------------------------------#
Analysis "protocol"
For all plots to do with pathways : python3 /Users/billiemeadowcroft/Dropbox/PHD_data/Dumbbells/dumbbellpathway.py output.xyz 1.0 
Which produces a file theta_d.dat in that directory
cp theta_d_${init}_${rat}_${seed}_${Eps}_${size}.dat
For pathway plot: 30x seeds with these theta files
#------------------------------------------------------------------------------#



#------------------------------------------------------------------------------#
For all plots to do with time to wrape (VS tensions) : python3 /Users/billiemeadowcroft/Dropbox/PHD_data/Dumbbells/wrapping.py output.xyz wrapping_${eps}_${R}_init_${init}_rat_${rat}_seed_${seed}.dat ${R}

Where R is the size of the dumbbell 
You need the Ovito module to run these. 
#------------------------------------------------------------------------------#