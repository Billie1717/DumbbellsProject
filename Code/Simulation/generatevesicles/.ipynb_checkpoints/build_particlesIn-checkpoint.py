from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
import math
import time
import random
# from scipy.spatial import KDTree
import sys, getopt
import os
# from numpy import cross, eye, dot
# from scipy.linalg import expm3, norm
# from scipy.optimize import minimize
# from mpl_toolkits.mplot3d import Axes3D

#'''
#This program generates randomly distributed spherical membrane-model particles in a box. It generates the Lammps read-in file.
#'''


#'''
#sys.argv [run_time] [init conc of solute] [ratio between concs- bigger than 1 = more on inside] [ratio between radius] [cores] [seed] [epsilon between cargo and membrane]
#'''


def cytoplasm(packing, cell_radius, cyto_radius, lbox, c_c0, Nmem, x_mem, y_mem, z_mem):
    print('Computing cytoplasm particles...')
    x_c = []
    y_c = []
    z_c = []
    x_in = []
    y_in = []
    z_in = []
    x_out = []
    y_out = []
    z_out = []
    nx = int(lbox/cyto_radius)                               # number of cyto particles per row
    ny = int((2*lbox)/(np.sqrt(3)*cyto_radius))
    nz = int((2*lbox)/((2 * np.sqrt(6)*cyto_radius)/3))
    # initialising a box of cytoplasm particles in centre of cell
    x_adjust = 0                                             # adjustments, making sure the cyto 'box' is in the middle of the cell
    y_adjust = 0 #11 * cyto_radius
    z_adjust = 0 #20 * cyto_radius
    eobx = -lbox + cyto_radius
    eoby = -lbox + np.sqrt(3) * cyto_radius/2
    eobz = -lbox + cyto_radius * np.sqrt(6)/3
    for k in range(nz):
        z_val = eobz + z_adjust + cyto_radius * (2 * np.sqrt(6) * k)/3
        for i in range (ny):
            if (i%2 == 0):
                for j in range (nx):
                    x_c.append(eobx + x_adjust + 2 * cyto_radius * j)
                    y_c.append(eoby + y_adjust + np.sqrt(3) * cyto_radius * i)
                    z_c.append(z_val)
            if (i%2 == 1):
                for j in range (nx):
                    x_c.append(eobx + x_adjust  + 2 * cyto_radius * j + cyto_radius)
                    y_c.append(eoby + y_adjust  + np.sqrt(3) * cyto_radius * i)
                    z_c.append(z_val)

    x_c = np.array(x_c)
    y_c = np.array(y_c)
    z_c = np.array(z_c)
    min_radius = cell_radius -  5 * cyto_radius               # making sure cyto doesn't touch membrane
    max_radius = cell_radius +  5 * cyto_radius               # making sure cyto doesn't touch membrane
    V_cell_in = (4/3)*np.pi*min_radius**3                     # volume of inside of cell
    V_cell_out = (4/3)*np.pi*max_radius**3                    # volume of inside of cell
    V_cyto = (4/3)*np.pi*cyto_radius**3                       # volume per cyto particle
    V_box = (2*lbox)**3
    num_in = int(np.rint(c_c0*packing*V_cell_in/V_cyto))
    num_out = int(np.rint(packing*(V_box-V_cell_out)/V_cyto))
    # RRR = np.sqrt(x_c**2 + y_c**2 + z_c**2)
    # max_nin = len(RRR[RRR<min_radius])
    # max_nout = len(RRR[RRR>max_radius])
    # max_cin = max_nin*V_cyto/V_cell_in
    # max_cout = max_nout*V_cyto/(V_box-V_cell_out)
    # print(max_cin,max_cout)
    m = 0
    n_i = 0
    n_o = 0
    av = np.zeros(len(x_c))
    idx = np.arange(len(av))
    while m < (num_in+num_out):
        if av.mean()==1:
            print('Maximum number of particles generated...')
            break
        # p = np.random.randint(len(x_c))
        p = random.choice(idx[av==0])
        if av[p] == 1:
            continue
        r = np.sqrt (x_c[p] **2 + y_c[p] **2 + z_c[p] **2)
        if r < min_radius and n_i < num_in:
            x_in.append(x_c[p])
            y_in.append(y_c[p])
            z_in.append(z_c[p])
            av[p] = 1
            n_i += 1
            m += 1
        elif r > max_radius and n_o < num_out:
            x_out.append(x_c[p])
            y_out.append(y_c[p])
            z_out.append(z_c[p])
            av[p] = 1
            n_o += 1
            m += 1
        else:
            av[p] = 1
        if m==np.rint(0.1*(num_in+num_out)):
            print('10%% of particles generated')
        elif m==np.rint(0.5*(num_in+num_out)):
            print('50%% of particles generated')
        elif m==np.rint(0.75*(num_in+num_out)):
            print('75%% of particles generated')
    actual_packing_in = len(x_in) * V_cyto /V_cell_in                      # packing of the cytoplasm in the cell
    actual_packing_out = len(x_out) * V_cyto /(V_box-V_cell_out)           # packing of the cytoplasm outside the cell
    print('Packing inside: %f\nPacking outside: %f'%(actual_packing_in,actual_packing_out))
    return x_in, y_in, z_in, x_out, y_out, z_out


def main():

    ################ Parameters ################
    # particle and concentration parameters

    c0 = float(sys.argv[2])             # initial concetration of solutes (homogeneous)
    c_c0 = float(sys.argv[3])           # ratio between concentrations
    r_r0 = c_c0**(1/float(3))           # ratio between radius
    radius = 0.5 						# particle radius
    Vpart = 4/3.0 * math.pi * radius**3 # particle volume
    Nparts = 0                          # initialize to 0, compute later...
    Nbonds = 0							# total # of bonds
    Nparticletypes = 4		# mem & inside & outside			# total # of particle types
    Nbondtypes = 0						# total # of bond types
    # membrane parameters
    rc_global = 2.6
    rc = 2.6
    rmin = 1.12
    mu = 3
    zeta = 4
    eps = 4.34
    sigma = 1.00
    theta0_11 = 0
    eps_mp = 2
    sigma_mp = sigma
    rc_mp = sigma * 2**(1/6.0)
    eps_pp = 2
    sigma_pp = sigma
    rc_pp = sigma * 2**(1/6.0)
    
    ################ cargo parameters ################
    sigma_c = float(sys.argv[8])
    eps_mc = int(sys.argv[7])
    sigma_mc = (sigma+sigma_c)/2
    rc_mc = (1.122462)*sigma_mc*1.2
    eps_pc = 2
    sigma_pc = (sigma_pp+sigma_c)/2
    rc_pc = (1.122462)*sigma_pc*1.2
    
    # simulation parameters
    seed = int(sys.argv[6])
    np.random.seed(seed)              					# random seed for permutations
    lbox = 100						                	# half box length
    Vbox = (2*lbox)**3                                  # box volume
    run_time = int(sys.argv[1])					        # total # of simulation time steps
    pbc = 1							                    # periodic boundary conditions
    dt = 0.01							                # time step
    thermo_interval = 2000                              # interval for measuring thermo quants.
    dump_interval = 2000                                # interval for saving coordinates
    # lammps submit file parameters
    cores = int(sys.argv[5])

    
    ################ Build initial configuration ################

    ################ Spherical Membrane ################
    Nmem = 4322
    workingdir = '/home/ucapbbm/Scratch/Dumbbells/scripts/'
    fname_in = workingdir +'RelaxedVesicalFiles/%d_relaxed_vesicle.xyz'%(Nmem) 			# read in positions of membrane particles
    f_in = open(fname_in,'rU')
    tmp = f_in.readline()
    tmp = f_in.readline()					#? two times necessary? YES - skip 2 lines
    #lbox = int(tmp.split()[1])
    lbox = int(tmp.split()[1])+12
    x_mem = []
    y_mem = []
    z_mem = []
    r_mem = []
    while True:
        tmp = f_in.readline()
        if tmp:
            x_mem.append(float(tmp.split()[1]))
            y_mem.append(float(tmp.split()[2]))
            z_mem.append(float(tmp.split()[3]))
            r_mem.append(np.sqrt(float(tmp.split()[1])**2+float(tmp.split()[2])**2+float(tmp.split()[3])**2))
        else:
            break

    cell_radius = np.mean(r_mem)				# radius of the cell (-1 such that the generated filament definitely is inside the cell)
    Nmem = len(x_mem)
    zMAX = max(z_mem)
    print(zMAX)
    xCM = np.sum(x_mem)/float(Nmem)
    yCM = np.sum(y_mem)/float(Nmem)
    zCM = np.sum(z_mem)/float(Nmem)
    print('Center of mass: [%.3f,%.3f,%.3f]'%(xCM,yCM,zCM))
    x_mem = x_mem - xCM
    y_mem = y_mem - yCM
    z_mem = z_mem - zCM
    xCM = np.sum(x_mem)/float(Nmem)
    yCM = np.sum(y_mem)/float(Nmem)
    zCM = np.sum(z_mem)/float(Nmem)
    print('Cell radius: %.3f'%(cell_radius))
    print('Center of mass: [%.3f,%.3f,%.3f]'%(xCM,yCM,zCM))

    ################ Outside and Inside particles ################
    cyto_radius = float(sys.argv[4])*1.12                                # radius of cytoplasm particles
    x_in, y_in, z_in, x_out, y_out, z_out = cytoplasm (c0, cell_radius, cyto_radius, lbox, c_c0, Nmem, x_mem, y_mem, z_mem)      # generates cytoplasm
    Nin = len(x_in)
    Nout = len(x_out)
    Ndumb=2
    Ntot = Nmem + Nin + Nout+Ndumb
    
    print('Total # of particles in the system: %d\n\tInside: %d\n\tOutside: %d\n\tMembrane: %d\n\tDumbbell: %d'%(Ntot,Nin,Nout,Nmem,Ndumb))

    ################  Write particles.in ################
    print('Writing particles.in\n')
    os.mkdir(workingdir + 'runs/SizeEps/init_{:.3f}_{:.3f}_{:d}_{:d}_{:.1f}'.format(c0,c_c0,seed,eps_mc,sigma_c))
    outfile=workingdir + 'runs/SizeEps/init_{:.3f}_{:.3f}_{:d}_{:d}_{:.1f}/'.format(c0,c_c0,seed,eps_mc,sigma_c) + 'particles.in'
    f=open(outfile,'w')
    f.write('LAMMPS data file generated by tc387 with a python script for hex membrane\n')
    f.write('\n')
    f.write(str(Ntot) + '   atoms\n')
    f.write(str(Nbonds) + '   bonds\n')
    f.write('\n')
    f.write(str(Nparticletypes)+' atom types\n')
    f.write(str(Nbonds)+' bond types\n')
    f.write('\n')
    f.write(str(-lbox)+' '+str(lbox)+' xlo xhi \n')
    f.write(str(-lbox)+' '+str(lbox)+' ylo yhi \n')
    f.write(str(-lbox)+' '+str(lbox)+' zlo zhi \n')
    f.write('\n')
    f.write('Masses\n')
    f.write('\n')
    for ii in range(Nparticletypes):
        f.write(str(ii+1) + ' 1\n')
    f.write('\n')
    f.write('Atoms # hybrid\n')
    f.write('\n')
    for ii in range(Nmem):
        r = float(np.sqrt(x_mem[ii]**2+y_mem[ii]**2+z_mem[ii]**2))
        f.write(str(ii+1) + '\t' + '1' + '\t' + str(x_mem[ii]) + '\t' + str(y_mem[ii]) + '\t' + str(z_mem[ii]) + '\t1\t1\t0\t' +   str(x_mem[ii]/r) + '\t' + str(y_mem[ii]/r) + '\t' + str(z_mem[ii]/r) + '\t' + str(ii+1) + '\n')
    for ii in range(Nin):
        f.write(str(1+Nmem+ii) + '\t' + '2' + '\t' + str(x_in[ii]) + '\t' + str(y_in[ii]) + '\t' + str(z_in[ii]) + '\t1\t1\t0\t0\t0\t0\t' + str(ii+1) + '\n')
    for ii in range(Nout):
        f.write(str(1+Nmem+Nin+ii) + '\t' + '3' + '\t' + str(x_out[ii]) + '\t' + str(y_out[ii]) + '\t' + str(z_out[ii]) + '\t1\t1\t0\t0\t0\t0\t' + str(ii+1) + '\n')
    for ii in range(Ndumb):
        f.write(str(Nout+Nmem+Nin+ii+1) + '\t' + '4' + '\t' + str(0) + '\t' + str(0) + '\t' + str(zMAX+ii*(sigma_c)+sigma_c*0.5+0.1) + '\t1\t1\t0\t0\t0\t0\t' + str(Nout+Nmem+Nin+1) + '\n')
    f.write('\n')
    f.close()
        

    ################ write in.local ################
    print('Writing in.local\n')
    outfile2 = workingdir + 'runs/SizeEps/init_{:.3f}_{:.3f}_{:d}_{:d}_{:.1f}/in.local'.format(c0,c_c0,seed,eps_mc,sigma_c)
    fo=open(outfile2,'w')
    fo.write('units\tlj\n')								#choice of reduced units
    fo.write('atom_style\thybrid sphere dipole molecular\n')				#choice of atom style (need spherical and dipole and molecular)
    fo.write('\n')
    fo.write('dimension\t3\n')								#working in 3D
    fo.write('boundary\tp p p\n')								#PBCs
    fo.write('processors\t* * 1\n')							#set the # of MPI processors in Lz dir. to  1	
    fo.write('\n')
    fo.write('read_data\t"particles.in"\n')						#read particle data from file
    fo.write('\n')
    fo.write('group\tmem\ttype 1\n')							#group membrane particles
    fo.write('group\tpar\ttype 2:3\n')							#group membrane particles
    fo.write('group\trest\ttype 1:3\n')							#group membrane particles
    fo.write('group\tcargo\ttype 4\n')							#group membrane particles
    fo.write('\n')
    fo.write('set group all mass 1.0\n')
    fo.write('\n')
    fo.write('################ membrane parameters ################\n')
    fo.write('variable\trc_global\tequal\t'+str(rc_global)+'\n')
    fo.write('variable\trc\tequal\t'+str(rc)+'\n')
    fo.write('variable\trmin\tequal\t'+str(rmin)+'\n')
    fo.write('variable\tmu\tequal\t'+str(mu)+'\n')
    fo.write('variable\tzeta\tequal\t'+str(zeta)+'\n')
    fo.write('variable\teps\tequal\t'+str(eps)+'\n')
    fo.write('variable\tsigma\tequal\t%.2f\n'%(sigma))
    fo.write('variable\ttheta0_11\tequal\t'+str(theta0_11)+'\n')
    fo.write('\n')
    fo.write('################ particle interaction parameters ################\n')
    fo.write('variable\teps_mp\tequal\t'+str(eps_mp)+'\n')
    fo.write('variable\tsigma_mp\tequal\t'+str(sigma_mp)+'\n')
    fo.write('variable\trc_mp\tequal\t'+str(rc_mp)+'\n')
    fo.write('variable\teps_pp\tequal\t'+str(eps_pp)+'\n')
    fo.write('variable\tsigma_pp\tequal\t'+str(sigma_pp)+'\n')
    fo.write('variable\trc_pp\tequal\t'+str(rc_pp)+'\n')
    fo.write('\n')
    fo.write('################ membrane and cargo ################\n')
    fo.write('variable\tsigma_c\tequal\t'+str(sigma_c)+'\n')
    fo.write('variable\teps_mc\tequal\t'+str(eps_mc)+'\n')
    fo.write('variable\tsigma_mc\tequal\t'+str(sigma_mc)+'\n')
    fo.write('variable\trc_mc\tequal\t'+str(rc_mc)+'\n')
    
    fo.write('################ particle and cargo ################\n')
    fo.write('variable\teps_pc\tequal\t'+str(eps_pc)+'\n')
    fo.write('variable\tsigma_pc\tequal\t'+str(sigma_pc)+'\n')
    fo.write('variable\trc_pc\tequal\t'+str(rc_pc)+'\n')
    fo.write('\n')
    fo.write('variable\tvatom\tequal\t4/3*PI*'+str(sigma)+'^3\n')
    fo.write('variable\tpatom atom -(c_stress[1]+c_stress[2]+c_stress[3])/v_vatom\n')
    fo.write('\n')
    fo.write('################ pair style ################\n')
    fo.write('# use hybrid overlay for pair_style membrane and lj/cut\n')
    fo.write('pair_style\thybrid/overlay\tmembrane ${rc_global} lj/cut ${rc_global}\n')
    fo.write('#pw359: Initialise LJ/expand to zero for all possible combinations\n')
    fo.write('pair_coeff\t*\t*\tlj/cut\t0\t0\n')
    fo.write('pair_coeff\t1\t1\tmembrane\t${eps}\t${sigma}\t${rmin}\t${rc}\t${zeta}\t${mu}\t${theta0_11}\n')
    fo.write('pair_coeff\t1\t2\tlj/cut\t${eps_mp}\t${sigma_mp}\t${rc_mp}\n')		            # volume exclusion membrane and inside parts
    fo.write('pair_coeff\t1\t3\tlj/cut\t${eps_mp}\t${sigma_mp}\t${rc_mp}\n')                    # volume exclusion membrane and outside parts
    fo.write('pair_coeff\t2\t2\tlj/cut\t${eps_pp}\t${sigma_pp}\t${rc_pp}\n')					# volume exclusion among cyto parts
    fo.write('pair_coeff\t2\t3\tlj/cut\t${eps_pp}\t${sigma_pp}\t${rc_pp}\n')                    # volume exclusion among cyto parts
    fo.write('pair_coeff\t3\t3\tlj/cut\t${eps_pp}\t${sigma_pp}\t${rc_pp}\n')                    # volume exclusion among cyto parts
    fo.write('pair_coeff\t1\t4\tlj/cut\t${eps_mc}\t${sigma_mc}\t${rc_mc}\n')					# volume exclusion among cyto parts
    fo.write('pair_coeff\t2\t4\tlj/cut\t${eps_pc}\t${sigma_pc}\t${rc_pc}\n')                    # volume exclusion among cyto parts
    fo.write('pair_coeff\t3\t4\tlj/cut\t${eps_pc}\t${sigma_pc}\t${rc_pc}\n')                    # volume exclusion among cyto parts
    fo.write('pair_modify\tshift\tyes\n')
    fo.write('\n')
    fo.write('################ neigh modify ################\n')
    fo.write('# Reduce the delay from default 10 to 2 to get rid of dangeours builds\n')
    fo.write('neigh_modify\tdelay 2\n')
    fo.write('neigh_modify\tpage 200000 one 20000\n')
    fo.write('\n')
    fo.write('comm_modify\tcutoff 50\n')
    fo.write('\n')
    fo.write('################ integrators ################\n')
    fo.write('fix\tfLANG\tall\tlangevin 1.0 1.0 1 ' +str(seed)+ ' zero yes omega yes\n')
    fo.write('fix\tfNVE\trest\tnve/sphere update dipole\n')
    fo.write('\n')
    fo.write('compute\tstress all stress/atom NULL\n')
    fo.write('compute\tcluster_mem mem cluster/atom 2.5\n')
    fo.write('compute\tcluster_cyto par cluster/atom 2.5\n')
    fo.write('\n')
    fo.write('################ dump ################\n')
    fo.write('dump\tcoords\tall\tcustom '+str(int(dump_interval/4))+' output.xyz id mol type x y z mux muy muz c_stress[1] c_stress[2] c_stress[3] c_stress[4] c_stress[5] c_stress[6] v_patom c_cluster_mem c_cluster_cyto\n')
    fo.write('dump_modify\tcoords sort id\n')
    fo.write('\n')
    fo.write('restart '+str(dump_interval)+'\tcontinue.dat\tcontinue.dat\n')
    fo.write('\n')
    fo.write('timestep\t'+str(dt)+'\n')
    fo.write('thermo\t'+str(thermo_interval)+'\n')
    fo.write('run\t2000\n')
    fo.write('fix dumbell cargo rigid/nve/small molecule\n')
    fo.write('################ run ################\n')
    fo.write('timestep\t'+str(dt)+'\n')
    fo.write('thermo\t'+str(thermo_interval)+'\n')
    fo.write('run\t'+str(run_time)+'\n')
    fo.close()
    
if __name__ == "__main__":
    main()
