# Preprocessing and postprocessing tools

This directory contains a collection of helpful preprocessing and
postprocessing utilities.

## Preprocessing scripts

### naluwind2amrwind.py
The [naluwind2amrwind.py](naluwind2amrwind.py) python script allows you to
convert a Nalu-wind input file to an amr-wind input file.  To the best of 
its ability, it will try to take the parameters from Nalu-wind and map them 
to amr-wind.

The basic usage is given by:  
```bash
usage: naluwind2amrwind.py [-h] [--outfile OUTFILE] yamlfile

Convert Nalu-Wind input file to AMR-wind input file

positional arguments:
  yamlfile

optional arguments:
  -h, --help         show this help message and exit
  --outfile OUTFILE  write output to this file
```

For instance, to convert the Nalu-wind input file `naluwind.yaml`, run
```bash
$ ./naluwind2amrwind.py naluwind.yaml

#---------------------------------------#
#            SIMULATION STOP            #
#---------------------------------------#
time.stop_time                 = 20000.0                        # Max (simulated) time to evolve
time.max_step                  = 40000                          # Max number of time steps

#---------------------------------------#
#         TIME STEP COMPUTATION         #
#---------------------------------------#
time.fixed_dt                  = 0.5                            # Use this constant dt if > 0
time.cfl                       = 0.95                           # CFL factor

#---------------------------------------#
#            INPUT AND OUTPUT           #
#---------------------------------------#
time.plot_interval             = 1000                           # Steps between plot files
time.checkpoint_interval       = 1000                           # Steps between checkpoint files
io.KE_int                      = 1                              
io.line_plot_int               = 1                              

[...snip...]
```
Use the `--outfile` option to write the output to a file instead.  Note 
that the `nalu_abl_mesh` and `nalu_preprocess` section should be present 
in the yaml file for it to correctly extract the geometry and the 
inversion layer properties.

## Postprocessing scripts

### postproamrwind.py

The [postproamrwind.py](postproamrwind.py) python provides a quick
method to time-average and extract data from the `line_plot.txt`
output.

The major functions in this code are:
- `loadfile(filename, t1, t2, headerinfo=[])`: Loads the line_plot.txt
  file given by `filename`, for all times between `t1` and `t2`.
- `timeaverage(dat, t1, t2, Ncells, Ncomp)`: Time averages the
  line_plot data given by `dat`, between times `t1` and `t2`.
- `avglineplot(filename, t1,t2, headerinfo=[])`: Combines `loadfile()`
  and `timeaverage()`, so it loads the line_plot.txt given by
  `filename`, and averages it between `t1` and `t2`.
- `extractvars(dat, colheaders, varnames)`: Returns the variable(s)
  matching the variables in varnames

A short example of how to use the code is shown below.  In this case,
we load the `line_plot.txt` file, average from 0 to 100 seconds, and
plot the averaged U velocity versus z.

```python
import postproamrwind
import matplotlib.pylot as plt
avgdat, headers = postproamrwind.avglineplot('line_plot.txt',0,100)
amrvars         = postproamrwind.extractvars(avgdat, headers, ['z','u_avg'])
plt.plot(amrvars['u_avg'], amrvars['z'])
plt.show()
```
