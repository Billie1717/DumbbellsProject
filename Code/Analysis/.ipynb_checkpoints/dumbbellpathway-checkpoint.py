import sys
import math
import numpy as np
from scipy.optimize import minimize 

def parsememb(infl): # infl = open("stresstensor.xyz", 'r') or infl = open("output.xyz", 'r')
    steps = []
    geos = []
    ts = 0
    atomnum = 0
    while 1: # processing one step
        current = []
        line = infl.readline()
        if line.startswith("ITEM: TIMESTEP"):
            nextline = infl.readline().split()
            ts = int(nextline[0])

            nextline = infl.readline()
            assert nextline.startswith("ITEM: NUMBER OF ATOMS")
            nextline = infl.readline().split()
            atomnum = int(nextline[0])

            nextline = infl.readline()
            assert nextline.startswith("ITEM: BOX BOUNDS")
            nextline = infl.readline()
            nextline = infl.readline()
            nextline = infl.readline()

            nextline = infl.readline()
            # if parsing "output.xyz"
            if nextline.startswith("ITEM: ATOMS id mol type x y z"):
# read all atom lines for this step
                for i in range(atomnum):
                    nextline = infl.readline().strip('\n').split()
                    if int(nextline[2]) == 4: # type 1 is membrane
                        current.append(list(map(float, nextline[3:6]))) # x, y, z
            else:
                print("Error! File cannot be parsed!")
                break
        else: # fail to readline, must be EOF
            break
        steps.append(ts)
        geos.append(current)
    # shift steps
    shifted_steps = [item - steps[0] for item in steps]
    return shifted_steps,geos

def dist2D(a,b):
    dx = a[0] - b[0]
    dy = a[1] - b[1]
    return math.sqrt(dx**2 + dy**2)

def dist3D(a,b):
    dx = a[0] - b[0]
    dy = a[1] - b[1]
    dz = a[2] - b[2]
    return math.sqrt(dx**2 + dy**2 + dz**2)

# single snapshot 
def pathway(x1,x2,y1,y2,z1,z2):
    lx = x1-x2
    ly = y1-y2
    lz = z1-z2
    xcom = (x1+x2)/2
    ycom = (y1+y2)/2
    zcom = (z1+z2)/2
    lDOTn = lx*xcom +ly*ycom + lz*zcom
    L = np.sqrt(lx**2 + ly**2 + lz**2)
    N = np.sqrt(xcom**2 + ycom**2 + zcom**2)
    theta = np.arccos(lDOTn/(L*N))
    theta_ = np.pi- min([abs(theta), abs(theta - np.pi)])
    return theta_, N


def rmse(lst):
    aver = sum(lst)/len(lst)
    resid2 = [(item - aver)**2 for item in lst]
    return np.sqrt(sum(resid2)/len(resid2))


# single frame
def xlimits(memb): #enter coords of 1 frame, takes out the tube ends 
#    print memb[:10]
    xvals = [item[0] for item in memb] 
    xvals_sorted = sorted(xvals)
    xs = sum(xvals_sorted[:50])/50
    xl = sum(xvals_sorted[-50:])/50
    return xs, xl

def trim(geo, xs, xl):
    geo_t = []
    for line in geo:
        if line[0] >= xs and line[0] <= xl:
            geo_t.append(line)
    return geo_t
    

# main
# zmin, zmax: either input parameter
# or zmin= zmin(memb)+20, zmax = zmax(memb)-5    
def main():
    nargs = len(sys.argv) - 1
    #if nargs != 3:
    #    sys.exit("Need 3 arguments!\nUsage: python radius.py output.xyz -20.0 -10.0 1.0\n")
    if nargs != 2:
        sys.exit("Example usage: python radius.py output.xyz 1.0\n")

    filename = sys.argv[1]
    #zmin = float(sys.argv[2])
    #zmax = float(sys.argv[3])
    bs = float(sys.argv[2])

    ofl = open("dumbbellcoords.dat", 'w')
    ofl2 = open("theta_d.dat",'w')
# step 1: parse membrane
    timesteps, geo = parsememb(open(filename, 'r'))
    #timesteps, geo = lp.parsemembpartial(open(filename, 'r'), zmin, zmax)
    for i in range(len(timesteps)):#range(len(timesteps-30)):
        global coords
        global dipoles
        coords = geo[i] 
        #print(coords)
        xmin, xmax = xlimits(coords) 
        #print("xmin and max:",xmin,xmax)
        #coords = trim(coords, xmin, xmax) #taking out the ends of the tube?
        x1 = coords[0][0]
        x2 = coords[1][0]
        y1 = coords[0][1]
        y2 = coords[1][1]
        z1 = coords[0][2]
        z2 = coords[1][2]
        theta = pathway(x1,x2,y1,y2,z1,z2)[0]
        d = pathway(x1,x2,y1,y2,z1,z2)[1]
# step 4: bin and calculate radius along z 
        #nw = tube_width(xmax, xmin, coords, bs)[0]
        #sw = sphere_diameter(xmax,xmin,coords,bs)
        #nwstd = tube_width(xmax, xmin, coords, bs)[1]# single snapshot
        ofl.write("%.2f %.2f %.2f %.2f %.2f %.2f %.2f\n" %(i, coords[0][0],coords[0][1],coords[0][2],coords[1][0],coords[1][1],coords[1][2]))
        ofl2.write("%.2f %.2f %.2f\n" %(i, theta,d))
        #ofl2.write("%.2f %.2f\n" %(i, sw))
        
if __name__ == "__main__":
    main()

     
