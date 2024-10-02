# assume cargo = type 3 and lipid head = type 1, lipid tail  = type 2
# line 58, line 59 need change if system change (eg, change filament length)
import sys
from ovito.io import *
from ovito.modifiers import *
from ovito.io import export_file
from ovito.data import *
import numpy as np


# output is a list of all particle-particle wrapN time series
def wrapN(pipeline, rcut): 
    output = []
    # period boundary condition unwrap
    pipeline.modifiers.append(UnwrapTrajectoriesModifier()) 

    # coordination calculation
    modifier = ComputePropertyModifier(cutoff_radius = rcut, expressions = 'NumNeighbors', output_property = 'Coordination')
    pipeline.modifiers.append(modifier)
    #data = pipeline.compute()

    # execute pipeline
    for frame in range(pipeline.source.num_frames):
        data = pipeline.compute(frame)
        wrapping = data.particles['Coordination'] # sum of ESCRT coordinations 
        output.append(wrapping)
    return output

# output is a list of escrt particle-particle wrapN time series


def main():
    pipeline = import_file(sys.argv[1], multiple_frames = True)
    rcut = sys.argv[3] + 1.0 # cut-off of interaction/coordination analysis
    wrapping = np.array(wrapN(pipeline, rcut)) # - np.array(wrapN_escrt(pipeline, rcut, Nescrt))
    print(wrapping)
    ofl = open(sys.argv[2], 'w')
    for i in range(len(wrapping)):
        ofl.write("%s %s\n" %(i, wrapping[i]))


if __name__ == "__main__":
    main()