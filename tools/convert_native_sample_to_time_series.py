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
sp_index = 0
s_name = "* will use first *"

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
    s_name = sys.argv[5]
if (narg > 6):
    sp_index = sys.argv[6]
if (narg > 7):
    print("WARNING: routine only takes 6 arguments maximum. Ignoring additional arguments")

print("\nParameters to be used:")
print("   run directory path: ", os.path.realpath(loc_dir))
print("   AMR-Wind path: ", os.path.realpath(AMR_WIND_PATH))
print("   post-processing directory name: ", pp_dir)
print("   top-level sampling label name: ", sampling_label)
print("   bottom-level sampling label name: ", s_name)
print("   particle index within sampler: ", sp_index)
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

# Separate by length, then order again (because 10000 would be next to 100000)
max_len = len(sdirs[0])
min_len = max_len
for n in range(len(sdirs)):
    max_len = max(max_len,len(sdirs[n]))
    min_len = min(min_len,len(sdirs[n]))

sdirs_init = list(sdirs)
if (max_len != min_len):
    cur_len = min_len
    idx = 0
    for inc_len in range(max_len - min_len + 1):
        cur_len += inc_len
        for n in range(len(sdirs_init)):
            if (len(sdirs_init[n]) == cur_len):
                sdirs[idx] = sdirs_init[n]
                idx += 1

n_dir = len(sdirs)
max_step = int(sdirs[n_dir-1][len(sampling_label):len(sdirs[n_dir-1])])

# Get stuff for header and array sizes using first available directory
n0 = 0
step = int(sdirs[n0][len(sampling_label):len(sdirs[n0])])
pt = pfile.load(step, label=sampling_label, root_dir = ".")
pt.parse_header()
pt.parse_info()
output = np.zeros((n_dir, 1 + 3 + pt.num_reals + pt.num_ints - 3))
v_str = "x y z"
for rv in range(pt.num_reals):
    v_str += " " + pt.real_var_names[rv]
for iv in range(pt.num_ints - 3):
    v_str += " " + pt.int_var_names[iv + 3]
header = "t " + v_str
samplers = pt.info["samplers"]
s_index = 0
if (s_name[0] != "*"):
    for s in range(len(samplers)):
        if (s_name == samplers[s]["label"]):
            s_index = s
else:
    s_name = samplers[0]["label"]

fname = sampling_label + "_" + s_name + "_" + str(sp_index) + "_time_series.txt"
print("Output file name: ", fname)

# Loop through files to get time series data
for n in range(n_dir):
    step = int(sdirs[n][len(sampling_label):len(sdirs[n])])
    print(str(step) + " out of " + str(max_step))
    pt = pfile.load(step, label=sampling_label, root_dir = ".")
    pt.parse_header()
    pt.parse_info()
    pt.load_binary_data()
    data = pt.df
    
    output[n, 0] = pt.info["time"]
    
    for p in range(pt.num_particles):
        if (data.set_id[p] == s_index and data.probe_id[p] == sp_index):
            output[n, 1] = data.xco[p]
            output[n, 2] = data.yco[p]
            output[n, 3] = data.zco[p]
            for rv in range(pt.num_reals):
                output[n, 4+rv] = data[pt.real_var_names[rv]][p]
            for iv in range(pt.num_ints - 3):
                output[n, 4+pt.num_reals+iv] = data[pt.int_var_names[iv+3]][p]
                
np.savetxt(fname, output[0:len(sdirs)],header=header)