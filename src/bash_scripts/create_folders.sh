#!/bin/bash

folder=`sed -n 1"p" name_folders_file.dat`
folder_event=`sed -n 2"p" name_folders_file.dat`
ls $folder > Magnitude_ev
NUM_LINES=$(wc -l "Magnitude_ev" | awk '{print $1}')
NUM_SCALING=$(wc -l "config_files/Parameters/classes_scaling.dat" | awk '{print $1}')

mkdir $folder_event
mkdir $folder_event/homogeneous_mu           
mkdir $folder_event/variable_mu   
for j in `seq 1 $NUM_LINES`; do
   Mw=`sed -n $j"p" Magnitude_ev`
   #echo $Mw
   Mw_string=$(echo $Mw | tr "." "_")
   mkdir $folder_event/homogeneous_mu/$Mw_string
   mkdir $folder_event/variable_mu/$Mw_string
   for i in `seq 1 $NUM_SCALING`; do
	 CL=`sed -n $i"p" config_files/Parameters/classes_scaling.dat`   
         mkdir $folder_event/homogeneous_mu/$Mw_string/$CL
         mkdir $folder_event/variable_mu/$Mw_string/$CL 
   done
done
