#!/usr/bin/env python
#
# Script to convert from a Nalu-wind YAML file to an AMR-wind input file
#
#
# usage: naluwind2amrwind.py [-h] [--outfile OUTFILE] yamlfile
#
# positional arguments:
#   yamlfile
#
# optional arguments:
#   -h, --help         show this help message and exit
#   --outfile OUTFILE  write output to this file


# In case anybody tries to use print() in python2
from __future__ import print_function

import sys
import math
import argparse

# Load the appropriate yaml reader
try:
    import ruamel.yaml as yaml
    useruemel=True
    yaml = yaml.YAML()
except:
    import yaml as yaml
    useruemel=False

# Get rid of yaml warnings (if present)
try: 
    yaml.warnings({'YAMLLoadWarning': False})
except:
    stuff=0

# Return the dictionary corresponding to key
#  listdic = list of dictionaries, indexed by key
def getdict(listdic, key):
    for dic in listdic:
        if key in dic: return dic[key]
    return []

# Return the dictionary which contains a given key
#  listdic = list of dictionaries, one of which might have key
def getdict2(listdic, key):
    for dic in listdic:
        if key in dic: return dic
    return []

# Return the dictionary which has keyname equal to va
def getdicwithname(listdic, val, keyname='name'):
    for dic in listdic:
        if dic[keyname] == val: return dic
    print("Cannot find entry with name: "+name)
    sys.exit(1) # Should not get here
    return 

# ----------
# Write a parameter to the amr input file
# ----------
#  param:    parameter name in AMR-wind file
#  data:     data to be written
#  outfile:  file handle (or sys.stdout)
#  isstring: if True, will wrap quotes around the string data
#  comment:  optional comment after data
#  commentonly: if True, don't write any data, just the comment
#  prefix:   Add prefix to the beginning of the line 
def writeAMRparam(param, data, outfile,
                  isstring=False, comment="", commentonly=False, 
                  prefix=''):
    if commentonly:
        print(comment, file=outfile)        
    elif data != None:  # Do not write anything if data is empty
        # -- construct the payload string --
        if isinstance(data, list):   datalist=data
        else:                        datalist=[data]
        if isstring:     datastring = ' '.join('\"'+str(x)+'\"' for x in datalist)
        else:            datastring = ' '.join(str(x) for x in datalist)
        # -- add comments --
        if len(comment)>0:           commentstring="# "+comment
        else:                        commentstring=""
        # -- Add the full string --
        fullstring = prefix+'%-30s = %-30s %s'%(param, datastring, commentstring)
        print(fullstring, file=outfile)
    return
    

### ========== Set up files and input/output args  ==================

# Load the command line arguments
parser = argparse.ArgumentParser(description='Convert Nalu-Wind input file to AMR-wind input file')
parser.add_argument('--outfile', default='', help="write output to this file")
parser.add_argument('yamlfile') #, nargs='+')
args=parser.parse_args()

# Load the yaml file
infile    = args.yamlfile
with open(infile) as fp:
    #yamldata=yaml.load(fp, Loader=yaml.FullLoader)
    yamldata=yaml.load(fp)

# Set up the output stream
#  Either write to outfile, or write to screen
if len(args.outfile)>0:
    outfile=open(args.outfile, 'w')
else:
    outfile=sys.stdout

### ========== Start reading nalu yaml file ==================
### --- sim start/stop/dt ----
timeinputyaml=yamldata['Time_Integrators'][0]['StandardTimeIntegrator']
max_step = timeinputyaml['termination_step_count']
fixed_dt = timeinputyaml['time_step']
steptype = timeinputyaml['time_stepping_type']

stop_time=fixed_dt*max_step

### --- output and restart ----
outputyaml          = yamldata['realms'][0]['output']
try:    plot_interval       = outputyaml['output_frequency']
except: plot_interval       = 0
try:    checkpoint_interval = yamldata['realms'][0]['restart']['restart_frequency']
except: checkpoint_interval = 0

### --- Set the mesh and grid sizes ----
if 'nalu_abl_mesh' in yamldata:
    # -- mesh extents --
    meshvertices     = yamldata['nalu_abl_mesh']['vertices']
    x0               = meshvertices[0]
    x1               = meshvertices[1]    
    # -- mesh dimensions --
    meshdimensions   = yamldata['nalu_abl_mesh']['mesh_dimensions']
    # Check periodicity
    is_periodic          = [1,1,0]   # Change this later
else:
    # Use the defaults
    x0               = [0,    0,    0]
    x1               = [1000, 1000, 1000]
    meshdimensions   = [48, 48, 48]
    # Check periodicity
    is_periodic      = [1,1,0]

### -- Set physics ---
soloptsyaml = yamldata['realms'][0]['solution_options']
gravity   = getdict(soloptsyaml['options'], 'user_constants')['gravity']
density   = getdict(soloptsyaml['options'], 'user_constants')['reference_density']
temp      = getdict(soloptsyaml['options'], 'user_constants')['reference_temperature']
sourceterms = getdict(soloptsyaml['options'], 'source_terms')['momentum']
kappa     = getdict(soloptsyaml['options'], 'turbulence_model_constants')['kappa']
smagCs    = getdict(soloptsyaml['options'], 'turbulence_model_constants')['cmuCs']
laminar_prandtl    = getdict(soloptsyaml['options'], 'laminar_prandtl')['enthalpy']
turbulent_prandtl  = getdict(soloptsyaml['options'], 'turbulent_prandtl')['enthalpy']
turbmodel = soloptsyaml['turbulence_model']

# coriolis parameters
latitude  = getdict(soloptsyaml['options'], 'user_constants')['latitude']
try:    east_vector  = getdict(soloptsyaml['options'], 'user_constants')['east_vector']
except: east_vector  = None
try:    north_vector = getdict(soloptsyaml['options'], 'user_constants')['north_vector']
except: north_vector = None
try:    secperrev    = 2*math.pi/getdict(soloptsyaml['options'], 'user_constants')['earth_angular_velocity']
except: secperrev    = None

# Turbulence model
if turbmodel != "smagorinsky":
    print("AMR-wind only supports smagorinsky model (so far), not "+turbmodel)
    sys.exit(1)
else:
    turbmodel = "Smagorinsky"

# Get surface roughness
z0      = getdict2(yamldata['realms'][0]['boundary_conditions'], 'wall_boundary_condition')['wall_user_data']['roughness_height']
zlotemp = getdict2(yamldata['realms'][0]['boundary_conditions'], 'wall_boundary_condition')['wall_user_data']['reference_temperature'] - temp
# Get upper temperature gradient
zhiTgrad=getdict2(yamldata['realms'][0]['boundary_conditions'], 'abltop_boundary_condition')['abltop_user_data']['normal_temperature_gradient']

# viscosity
viscosityyaml = getdicwithname(yamldata['realms'][0]['material_properties']['specifications'], 'viscosity')
viscosity     = viscosityyaml['value']

# Get forcing terms
forcingterms =  getdict(soloptsyaml['options'], 'source_terms')['momentum']
# Remap the forcing terms to AMR-wind forcing terms
forcingdictmap = { 'buoyancy_boussinesq':'BoussinesqBuoyancy', 
                   'EarthCoriolis':'CoriolisForcing',
                   'abl_forcing':'ABLForcing'
               }
forcingterms = [forcingdictmap[x] for x in forcingterms]
#forcingterms = ["BoussinesqBuoyancy","CoriolisForcing","ABLForcing"]

# ABL forcing
ablyaml = yamldata['realms'][0]['abl_forcing']['momentum']
abl_forcing_height = ablyaml['heights']
abl_velx           = ablyaml['velocity_x'][0][1]
abl_vely           = ablyaml['velocity_y'][0][1]
abl_velz           = ablyaml['velocity_z'][0][1]
abl_velocity       = [abl_velx, abl_vely, abl_velz]

# temperature of inversion layer
if 'nalu_preprocess' in yamldata:
    tempheights = yamldata['nalu_preprocess']['init_abl_fields']['temperature']['heights']
    tempvals    = yamldata['nalu_preprocess']['init_abl_fields']['temperature']['values']
else:
    tempheights = [0, 0.5*(x1[2]-x0[2])+x0[2], x1[2]]
    tempvals    = [temp, temp, temp]
    
### --- set other defaults ----
# Some stuff that wouldn't be specified in nalu
cfl               = 0.95

# AMR defaults
AMRdefaults=[
    ['io.KE_int',        1, ''],
    ['io.line_plot_int', 1, ''],
]

# Physics defaults
physicsdefaults = [
    ['incflo.use_godunov', 1,    ''],
    ['incflo.physics',    'ABL', ''],
]

# verbose defaults
verbosedefaults = [
    ['incflo.verbose',    3, 'incflo.level'],
]

# tolerances and debug defaults
debugdefaults = [
    ['amrex.fpe_trap_invalid',  0, 'Trap NaNs'],          
    ['diffusion.mg_verbose', 0, ''],
    ['diffusion.mg_cg_verbose',  0, ''],
    ['diffusion.mg_rtol',  1.0e-6, ''],
    ['diffusion.mg_atol',  1.0e-8, ''],
    ['mac_proj.mg_rtol',  1.0e-6, ''],
    ['mac_proj.mg_atol',  1.0e-8, ''], 
    ['nodal_proj.mg_verbose',  0, ''],
    ['nodal_proj.mg_rtol',  1.0e-6, ''],
    ['nodal_proj.mg_atol',  1.0e-8, ''],
]

### ========== Start writing out AMR inputfile ==================
stopheader="""
#---------------------------------------#
#            SIMULATION STOP            #
#---------------------------------------#
"""
writeAMRparam([], [], outfile,   comment=stopheader[:-1], commentonly=True)
writeAMRparam('time.stop_time',  stop_time, outfile, comment="Max (simulated) time to evolve")
writeAMRparam('time.max_step',   max_step, outfile, comment="Max number of time steps")

timeheader="""
#---------------------------------------#
#         TIME STEP COMPUTATION         #
#---------------------------------------#
"""
writeAMRparam([], [], outfile,   comment=timeheader[:-1], commentonly=True)
writeAMRparam('time.fixed_dt',   fixed_dt, outfile, comment="Use this constant dt if > 0")
writeAMRparam('time.cfl',        cfl,      outfile, comment="CFL factor")

ioheader="""
#---------------------------------------#
#            INPUT AND OUTPUT           #
#---------------------------------------#
"""
writeAMRparam([], [], outfile,   comment=ioheader[:-1], commentonly=True)
writeAMRparam('time.plot_interval', plot_interval, outfile, comment="Steps between plot files")
writeAMRparam('time.checkpoint_interval', checkpoint_interval, outfile, comment="Steps between checkpoint files")
# Other defaults
for default in AMRdefaults:
    writeAMRparam(default[0], default[1], outfile, comment=default[2])

physicsheader="""
#---------------------------------------#
#               PHYSICS                 #
#---------------------------------------#
"""
writeAMRparam([], [], outfile,   comment=physicsheader[:-1],  commentonly=True)
writeAMRparam('incflo.gravity',  gravity, outfile, comment="Gravitational force (3D)")
writeAMRparam('incflo.density',  density, outfile, comment="Reference density")
writeAMRparam('transport.viscosity',         viscosity,         outfile)
writeAMRparam('transport.laminar_prandtl',   laminar_prandtl,   outfile)
writeAMRparam('transport.turbulent_prandtl', turbulent_prandtl, outfile)
writeAMRparam('turbulence.model',            turbmodel,         outfile)
writeAMRparam('Smagorinsky_coeffs.Cs',       smagCs,            outfile)

writeAMRparam([], [], outfile,   comment="\n# ABL forcing", commentonly=True)
writeAMRparam('ICNS.source_terms', forcingterms,                outfile)
writeAMRparam('BoussinesqBuoyancy.reference_temperature', temp, outfile)
writeAMRparam('ABLForcing.abl_forcing_height', abl_forcing_height, outfile)
writeAMRparam('incflo.velocity',               abl_velocity,    outfile)

writeAMRparam([], [], outfile,   comment="", commentonly=True)
writeAMRparam('ABL.temperature_heights', tempheights[1:], outfile)
writeAMRparam('ABL.temperature_values',  tempvals[1:],    outfile)
writeAMRparam('ABL.kappa',               kappa,           outfile)
writeAMRparam('ABL.surface_roughness_z0', z0,             outfile)

writeAMRparam([], [], outfile,   comment="# Coriolis forcing", commentonly=True)
writeAMRparam('CoriolisForcing.latitude',     latitude,          outfile)
writeAMRparam('CoriolisForcing.east_vector',  east_vector,       outfile)
writeAMRparam('CoriolisForcing.north_vector', north_vector,      outfile)
writeAMRparam('CoriolisForcing.rotational_time_period', secperrev,   outfile)


for default in physicsdefaults:
    writeAMRparam(default[0], default[1], outfile, comment=default[2])

amrheader="""
#---------------------------------------#
#        ADAPTIVE MESH REFINEMENT       #
#---------------------------------------#
"""
writeAMRparam([], [], outfile,   comment=amrheader[:-1],  commentonly=True)
writeAMRparam('amr.n_cell',      meshdimensions, outfile, comment="Grid cells at coarsest AMRlevel");
writeAMRparam('amr.max_level',   0,              outfile, comment="Max AMR level in hierarchy");

geomheader="""
#---------------------------------------#
#              GEOMETRY                 #
#---------------------------------------#
"""
writeAMRparam([], [], outfile,   comment=geomheader[:-1], commentonly=True)
writeAMRparam('geometry.prob_lo',      x0,             outfile, comment="Lo corner coordinates")
writeAMRparam('geometry.prob_hi',      x1,             outfile, comment="Hi corner coordinates")
writeAMRparam('geometry.is_periodic',  is_periodic,    outfile, comment="Periodicity x y z (0/1)")
writeAMRparam([], [], outfile,   comment="\n# Boundary conditions", commentonly=True)
writeAMRparam('zlo.type',             'wall_model',    outfile, isstring=True)
writeAMRparam('zlo.temperature',       zlotemp,        outfile)
writeAMRparam('zhi.type',             'slip_wall',     outfile, isstring=True)
writeAMRparam('zhi.temperature',       zhiTgrad,       outfile)

verbheader="""
#---------------------------------------#
#              VERBOSITY                #
#---------------------------------------#
"""
writeAMRparam([], [], outfile,   comment=verbheader[:-1], commentonly=True)
for default in verbosedefaults:
    writeAMRparam(default[0], default[1], outfile, comment=default[2])

debugheader="""
#---------------------------------------#
#              DEBUGGING                #
#---------------------------------------#
"""
writeAMRparam([], [], outfile,   comment=debugheader[:-1], commentonly=True)
writeAMRparam([], [], outfile,   comment="##possible debugging parameters", commentonly=True)
for default in debugdefaults:  # write defaults (as comments)
    writeAMRparam(default[0], default[1], outfile, comment=default[2], prefix='# ')

### ========== Finished writing out AMR inputfile ==================
outfile.close()
