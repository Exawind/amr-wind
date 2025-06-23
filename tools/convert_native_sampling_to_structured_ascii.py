#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
import os
AMR_WIND_PATH = 'amr-wind'
loc_dir = "."
sampling_label = "sampling"
pp_dir = "post_processing"

########## Input arguments from command line ##############
narg = len(sys.argv)
print("Total input arguments passed:", narg-1)

print("\nArguments passed:", end = " ")
for i in range(1, narg):
    print(sys.argv[i], end = " ")

if (narg > 1):
    loc_dir = sys.argv[1]
if (narg > 2):
    AMR_WIND_PATH = sys.argv[2]
if (narg > 3):
    pp_dir = sys.argv[3]
if (narg > 4):
    sampling_label = sys.argv[4]
if (narg > 5):
    print("WARNING: routine only takes 4 arguments maximum. Ignoring additional arguments")

print("\nParameters to be used:")
print("   run directory path: ", os.path.realpath(loc_dir))
print("   AMR-Wind path: ", os.path.realpath(AMR_WIND_PATH))
print("   post-processing directory name: ", pp_dir)
print("   top-level sampling label name: ", sampling_label)
###########################################################

sys.path.append(AMR_WIND_PATH+'/tools/')
import amrex_particle
from amrex_particle import AmrexParticleFile

pp_dir = loc_dir + "/" + pp_dir
pfile = AmrexParticleFile(loc_dir)

os.chdir(pp_dir)
contents = glob.glob(sampling_label + "*")
directories = [path for path in contents if os.path.isdir(path)]
sdirs = sorted(directories)

max_dir = len(sdirs)-1
max_step = int(sdirs[max_dir][len(sampling_label):len(sdirs[max_dir])])

for n in range(len(sdirs)):
    step = int(sdirs[n][len(sampling_label):len(sdirs[n])])
    print(str(step) + " out of " + str(max_step))
    pt = pfile.load(step, label=sampling_label, root_dir = ".")
    pt.parse_header()
    pt.parse_info()
    pt.load_binary_data()
    data = pt.df
    samplers = pt.info["samplers"]
    
    t_str = str(pt.info["time"])
    v_str = "x y z"
    
    for rv in range(pt.num_reals):
        v_str += " " + pt.real_var_names[rv]
    for iv in range(pt.num_ints - 3):
        v_str += " " + pt.int_var_names[iv + 3]

    header = "t = " + t_str + "\n" + v_str
        
    # Have to go with max possible size, subsampler size not stored
    output = np.zeros((pt.num_particles, 3 + pt.num_reals + pt.num_ints - 3))
    
    for s in range(len(samplers)):
        np_s = 0
        for p in range(pt.num_particles):
            if data.set_id[p] == s:
                np_s += 1
                pid = data.probe_id[p]
                output[pid, 0] = data.xco[p]
                output[pid, 1] = data.yco[p]
                output[pid, 2] = data.zco[p]
                for rv in range(pt.num_reals):
                    output[pid, 3+rv] = data[pt.real_var_names[rv]][p]
                for iv in range(pt.num_ints - 3):
                    output[pid, 3+pt.num_reals+iv] = data[pt.int_var_names[iv+3]][p]
                

        fname = sdirs[n] + "_" + samplers[s]["label"]+ ".txt"
        np.savetxt(fname,output[0:np_s,:],header=header)
