#!/bin/bash

event=`sed -n 1"p" name_folders_file.dat`
event_out=`sed -n 2"p" name_folders_file.dat`
zone_code=`sed -n 3"p" name_folders_file.dat`
numb_stoch=`sed -n 4"p" name_folders_file.dat`
variable_mu=`sed -n 5"p" name_folders_file.dat`
ls $event > Magnitude_ev
NUM_MAGN=$(wc -l "Magnitude_ev" | awk '{print $1}')
folder_in=$(pwd | awk '{print $1}')
NUM_SCALING=$(wc -l "config_files/Parameters/classes_scaling.dat" | awk '{print $1}')

for j in `seq 1 $NUM_MAGN`; do
#for j in `seq 1  $NUM_LINES`; do
  #echo $j
  Mw=`sed -n $j"p" Magnitude_ev`
  Mw=$(echo $Mw | tr "_" ".")
  Mw_string=$(echo $Mw | tr "." "_")
#DOM: INTRODURRE UN CICLO SULLE 4 CLASSI O SI LANCIANO 4 SCRIPT PER OGNI CLASSE?   
  for i in `seq 1 $NUM_SCALING`; do
   CL=`sed -n $i"p" config_files/Parameters/classes_scaling.dat`   #!!!!!! CHECK CLASS
   rigidity=variable_mu  #!!!!!!!! CHECK RIGIDITY
   echo magnitude=$Mw > input_magnitude
   folder=$folder_in/$event/$Mw_string/$CL  #!!!!!!CHECK AREA AND CLASS
   echo $folder_in"/config_files/Matrix_distances/"$zone_code"_matrix_distance.bin" > matrix_string.txt
   #echo $folder
   folder_seed=$folder_in/$event_out/homogeneous_mu/$Mw_string/$CL   #!!!!!!!!CHECK CLASS
   folder_out=$folder_in/$event_out/$rigidity/$Mw_string/$CL   #!!!!!!!!CHECK CLASS
   echo $folder_out
   cp bin/k223d.x $folder_out
   cp input_magnitude param_zone.dat matrix_string.txt $folder_out
   cd config_files
   cp Mesh/$zone_code*mesh* Rigidity/mu_$zone_code.dat $folder_out
   cp Parameters/param_var.dat $folder_out/param.dat
   cd ..
   cp $folder/QuakeArea*.dat $folder_out
   cp $folder/Slip_PDF*.dat $folder_out
   cp $folder/mu_Slip_aux*.dat $folder_out
   #cp $folder_seed/Seed_00* $folder_out
   #cp $folder_seed/Seed_01* $folder_out
   #cp $folder_seed/Seed_02* $folder_out
   #cp $folder_seed/Seed_03* $folder_out
   #cp $folder_seed/Seed_04* $folder_out
   #cp $folder_seed/Seed_05* $folder_out
   ls $folder_out/QuakeArea*.dat > index_scenario
   NUM_SCENARIO=$(wc -l "index_scenario" | awk '{print $1}')
   echo $NUM_SCENARIO $numb_stoch >> index_file.dat
   for l in `seq 1 $NUM_SCENARIO`; do
       file_event=`sed -n $l"p" index_scenario`
       #echo $file_event
       eventid=${file_event: -9:5}
       #echo $eventid
#DOM al primo run creo i seed per il caso omogeneo e poi per i casi successivi li vado a leggere nel folder2?     
     #NUM_STOCHASTIC=$(wc -l "list_seed" | awk '{print $1}')
       for i in `seq 1 $numb_stoch`; do
       INDEX=$i
       #echo $INDEX
       if [ $INDEX -lt 10 ]
           then
              STRING4FILE=00$INDEX
              STRING4FILE=${STRING4FILE//' '/''}
           elif [ $INDEX -lt 100 ]
           then
              STRING4FILE=0$INDEX
              STRING4FILE=${STRING4FILE//' '/''}
           else
              STRING4FILE=$INDEX
              STRING4FILE=${STRING4FILE//' '/''}
        fi
        if [ $INDEX -le 2 ]
           then
              numb_gauss=1
           elif [  $INDEX -gt 2 ] && [ $INDEX -le 4 ]
           then
              numb_gauss=2
           else
              numb_gauss=3
        fi
        indexd=$eventid 
        indexd+="_"
        indexd+=$STRING4FILE
        echo $indexd $numb_gauss >> index_file.dat
     done
    done
    mv index_file.dat $folder_out
    cd $folder_out
    ./k223d.x input=param.dat > output_file.txt
    rm QuakeArea*.dat mu_Slip_aux*.dat Slip_PDF*.dat
    cd ../../../../
#mkdir /nas/cat2/scala/italian_map/CalabrianArc/variable_mu/$Mw_string
#mkdir /nas/cat2/scala/italian_map/CalabrianArc/variable_mu/$Mw_string/$CL
#mv $folder_out/*  /nas/cat2/scala/italian_map/$folder_out 
 done
done
if [ $variable_mu -eq 1 ]
	   then
              mkdir input
	      mkdir output
	      mv $event_out output
	      mv $event input
	      rm Magnitude_ev *.txt index_scenario input_magnitude *.dat
fi
