#!/bin/sh

export LD_LIBRARY_PATH=/Users/osterhof/.root/lib

export config_geometry_s1=5730000
export config_geometry_s2=5730000
export config_geometry_theta=0.0109
export config_illumination_center=2.0
export config_illumination_factor=1.0
export config_illumination_file=""
export config_illumination_intens=1.0
export config_illumination_method="self"
export config_illumination_number=1
export config_illumination_shm=""
export config_illumination_width=0.1
export config_layers_material1="W"
export config_layers_material2="B4C"
export config_layers_number=60
export config_layers_ratio=0.5
export config_mirror_length=99120
export config_physic_wavelength=0.00008265
export config_shm_fieldout=56156
export config_shm_plotdata=56155
export config_simulation_gridpointss=4720000
export config_simulation_gridpointst=1200
export config_xocd_hostname="k-raum.org"
export config_xocd_port=42001



export config_vertical_coherence=9
export config_horiz_coherence=31
export config_distance_coherence=40000000
export config_step_coherence=1
export config_divergence_correction=0

export config_Fname_MLtopo="Surf_HH.dat"
export config_interpolation_rate=1
export config_detector_NxCAM=2048
export config_detector_NyCAM=2048
export config_detector_SxCAM=0.732421875
export config_detector_SyCAM=0.732421875
export config_FWHMpsf=1.5
export config_detector_DistanceMIN=0
export config_detector_DistanceMAX=990000
export config_distance_step=3

export config_SimulationFile="F-test2RM2.dat"
export config_TTfieldfile="WVF-test2RM2.dat"
export config_TTfeaturefile="TP-test2RM2.dat"


make -q || make || exit && echo

./takagi-taupin -cp $@ --flat

#python ./Propag/Pro_Fresnel.py "field_out.dat" "To_propagate.dat"
