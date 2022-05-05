# -*- coding: utf-8 -*-
""" 
Calculate relevant quantities from statistics of precursor ABL-LES
for inflow/outflow simulation

Developed by Ganesh Vijayakumar and modified by Alex Rybchuk 
"""

import sys
from pathlib import Path
import argparse, textwrap
import numpy as np
import netCDF4 as nc
import os
from shutil import move
from tempfile import mkstemp

class ABLInflowOutflowPrep:
    """Prep inputs for Inflow/Outflow simulation"""

    def __init__(self, stats_file, t_start, t_end, input_file):
        """
        Args:
            stats_file (path): Path to the file containing statistics during
                               time of interest
            t_start (float): Start time of inflow/outflow simulation
            t_end (float): End time of inflow/outflow simulation
            input_file (path): Path to the input file for the inflow/outflow
                                 simulation
        """
        self.stats_file = stats_file
        self.dset = nc.Dataset(Path(stats_file))
        self.t_start = t_start
        self.t_end = t_end
        if input_file != "":
            self.input_file = Path(input_file)

    def calc_runtime_stats(self):
        """Examine statistics file and calculate relevant quantities """

        time = self.dset['time'][:]
        t_filter = (time > self.t_start) & (time < self.t_end)

        if (not t_filter.any()):
            dset_tstart = time[0]
            dset_tend = time[-1]
            sys.exit(textwrap.dedent(f"""
            Statistics file time window does not cover desired run time of inflow/outflow simulation

            Desired start time   = {str(self.t_start)}
            Desired end time     = {str(self.t_end)}

            Time window in statistics file = {str(dset_tstart)} - {str(dset_tend)}
            """))

        self.abl_force_x = np.average(self.dset['abl_forcing_x'][t_filter])
        self.abl_force_y = np.average(self.dset['abl_forcing_y'][t_filter])

        mean_profiles = self.dset['mean_profiles']
        self.hvelmag = np.average(mean_profiles['hvelmag'][t_filter, 0])
        self.u_mean = np.average(mean_profiles['u'][t_filter, 0])
        self.v_mean = np.average(mean_profiles['v'][t_filter, 0])
        self.theta_mean = np.average(mean_profiles['theta'][t_filter, 0])

        h = self.dset['mean_profiles']['h'][:]
        avg_theta = np.average(mean_profiles['theta'][t_filter, :], 0)
        with open('avg_theta.dat','w') as f:
            f.write('{} \n'.format(avg_theta.size))
            for i,t in enumerate(avg_theta):
                f.write('{} {} \n'.format(h[i], t))
            
        print(textwrap.dedent(f"""                 
        Statistics file      = {self.stats_file}
        Desired start time   = {str(self.t_start)} 
        Desired end time     = {str(self.t_end)}   
        abl_forcing_x        = {str(self.abl_force_x)}  
        abl_forcing_y        = {str(self.abl_force_y)}  
        hvelmag              = {str(self.hvelmag)}      
        u_mean               = {str(self.u_mean)}       
        v_mean               = {str(self.v_mean)}       
        theta_mean           = {str(self.theta_mean)}   
        """))
        
    def check_input_file(self):
        """
        Does the input file have any lines that will be duplicated?
        """
        
        ## Read file
        with open(self.input_file) as f:
            lines = f.readlines()
        
        ## Parse arguments in file
        for line in lines:                
            # Warn if duplicate lines will exist
            warning_list = ['ABL.wall_shear_stress_type',
                            'ABL.inflow_outflow_mode',
                            'ABL.wf_velocity',
                            'ABL.wf_mag',
                            'ABL.wf_theta',
                            'BoussinesqBuoyancy.read_temperature_profile',
                            'BoussinesqBuoyancy.tprofile_filename',
                            'BodyForce.magnitude']
            for warn in warning_list:
                if warn in line:
                    print(textwrap.dedent(f'Warning: "{warn}" already present'))
                    
        ## Double check if 'ABLMeanBoussinesq' is present
        mean_bouss_flag = False
        for line in lines:
            if "ABLMeanBoussinesq" in line:
                mean_bouss_flag = True
        if not mean_bouss_flag:
            print(textwrap.dedent(f'Info: ABLMeanBoussinesq not in ICNS.source_terms'))

        
    def update_input_file(self):
        """
        Add statistical info to the outflow file 
        """
        print(textwrap.dedent(str(self.abl_force_x)))
        
        intro_line = "#--------- Additions by calc_inflow_stats.py ---------#\n"
        wall_shear_line = 'ABL.wall_shear_stress_type = "local"\n'
        wall_io_line = "ABL.inflow_outflow_mode = true\n"
        wf_vel_line = f"ABL.wf_velocity = {self.u_mean} {self.v_mean}\n"
        wf_velmag_line = f"ABL.wf_mag = {self.hvelmag}\n"
        wf_theta_line = f"ABL.wf_theta = {self.theta_mean}\n"
        bodyforce_mag_line = f"BodyForce.magnitude = {self.abl_force_x} {self.abl_force_y} 0.0\n"
        bouss_read_line = "BoussinesqBuoyancy.read_temperature_profile = true\n"
        bouss_profile_line = "BoussinesqBuoyancy.tprofile_filename = avg_theta.dat\n"
        outro_line = "#-----------------------------------------------------#\n"
        
        ### Modify input file
        written_flag = False  # Used to deal with duplicate "ICNS.source_terms"
        fh, abs_path = mkstemp()
        with os.fdopen(fh,'w') as new_file:
            with open(self.input_file) as old_file:
                for i, line in enumerate(old_file):
                    # Add modifications near the ICNS.source_terms line
                    if ('ICNS.source_terms' not in line) or written_flag:  # Don't change anything
                        new_file.write(line)
                    else:  # Add in new lines
                        new_file.write(line)                            
                        new_file.write(intro_line)
                        new_file.write(wall_shear_line)
                        new_file.write(wall_io_line)
                        new_file.write(wf_vel_line)
                        new_file.write(wf_velmag_line)
                        new_file.write(wf_theta_line)
                        new_file.write(bodyforce_mag_line)
                        new_file.write(bouss_read_line)
                        new_file.write(bouss_profile_line)
                        new_file.write(outro_line)
                        
                        written_flag = True
        # Remove original file
        os.remove(self.input_file)
        # Move new file
        move(abs_path, self.input_file)

        
def main():
    """Run program"""
    parser = argparse.ArgumentParser(
        description="Prep ABL simulations for inflow/outflow turbine simulations")
    parser.add_argument(
        '-sf', '--statsfile', type=str, required=True,
        help="Name of the statistics netCDF file from precursor run")
    parser.add_argument(
        '-ts', '--t-start', type=float, required=True,
        help="Start time")
    parser.add_argument(
        '-te', '--t-end', type=float, required=True,
        help="End time")
    parser.add_argument(
        '-if', '--inputfile', required=False,
        help="[Optional] Name of the input file for the inflow/outflow simulation",
        default="")

    args = parser.parse_args()

    abl_if = ABLInflowOutflowPrep(
        args.statsfile, args.t_start, args.t_end, args.inputfile)
    abl_if.calc_runtime_stats()
    if args.inputfile != "":
        abl_if.check_input_file()
        abl_if.update_input_file()

if __name__ == "__main__":
    main()
