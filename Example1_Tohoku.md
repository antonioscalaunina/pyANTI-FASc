# EXAMPLE 1 Tohoku earthquake

In this brief guide, a practical example for running ANTI-FASc is shown. 
In this example, a set of slip distributions based on location and magnitude of 2011, March the 11th, Tohoku earthquake with a magnitude Mw=9.1, is generated.


# 1 - Mesh

pyANTI-FASc makes available an ensemble of predefined mesh discretizations for most of the seismogenic portion of the subducting plates worldwide. These meshes have been defined starting from the geometry defined in the framework of the project Slab 2.0 and available at this [webpage](https://www.sciencebase.gov/catalog/item/5aa1b00ee4b0b1c392e86467). They are discretized using a pretty homogeneous minimum nodes inter-distance size ranging from about 10 to 15 km. All the available meshes can be found [here](https://github.com/antonioscalaunina/pyANTI-FASc/tree/main/utils/sz_slabs) and the Figure 1 shows them on the map.

![Map of the subducting plate meshes available in the current version of pyANTI-FASc](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/utils/sz_slabs/map_of_slabs.png)
*Figure 1 - Map of the subducting plate meshes available in the current version of pyANTI-FASc*


# 2 - Input files

## 2.1 input.json

This file [input.json](https://github.com/antonioscalaunina/pyANTI-FASc/tree/main/config_files/Parameters/input.json) contained in the [config_files/Parameters](https://github.com/antonioscalaunina/pyANTI-FASc/tree/main/config_files/Parameters) folder manages the main configuration parameters for the run of ANTI-FASc. The file actually used for the run must be always named input.json. The default available file is set to run this example, but it can be easily modified to run other possible cases.

Here below the important settings to be managed by the user are shown. Look carefully to the comments besides the parameters. The parameters not shown in this example might be left unmodified. Their use and functionality will be fixed in next releases and better described in [Wiki Documentation](https://github.com/antonioscalaunina/pyANTI-FASc/wiki) currently under construction.

	{"zone_name": "kurilsjapan",      # Name of the precomputed mesh to be used. See the slab database made available within the repository and use the Mesh Folder for the slab you want to select
	"Merc_zone": 54,                  # Mercator zone for the selected slab. See the slab database and use the correct Mercator zone 
        "acronym": "KuJ",		  # 3 digit acronym that is used for that slab. It can be arbitrarily chosen by the user (but must have 3 digits!). You might find suggestions into the slab database 

The list of the slab database can be found [here](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/utils/sz_slabs/slabs_database) 
       	
	"Event": {
	"Name": "Tohoku_test",                    # Name of the test. It is used to define the output folder name
	"Hypo_LonLat" : [142.369, 38.322],        # Epicenter of the event
	"Magnitude" : 9.0                         # Magnitude of the event
	},
	"Configure": {
	"application": "PTF",                     # This application restricts the computed scenarios to a range of magnitude and location around predefined values
	"shape": "Rectangle",                     # This choice allows to compute scenarios with aspect ratio L/W preserved as prescribed by the selected scaling law. The other possible choice is "Circle". More details soon in the Wiki Documentation
	"numb_stoch": 5,                          # Number of stochastic slip for each rupture areas
	"variable_mu": 1,                         # 1 means that also the distributions with variable rigidity will be computed. 0 for computing only the case with homogeneous rigidity
	"coupling_shallow_limit":2.5,             # Shallow limit for the area where the seismic coupling is expected to decrease
	"coupling_deep_limit":40.0,		  # Deep limit for the area where the seismic coupling is expected to decrease


 
	"Magnitude_lb": 0.15,                    # Magnitude in a range [Mw-0.15 Mw+0.15] will be accounted, used only for "application": "PTF"
	"Magnitude_ub": 0.15,
        "hypo_baryc_distance": 1.0,             # Rupture barycenters at less than 1 Length form hypocenter will be used to define areas and slip distributions, used only for "application": PTF. The Lengthis inferred from scaling law for each Magnitude bin.
	"minimum_bnd_distance": 0.25,           # This option (as well as the next one) is used to limit the number of rupture areas dependending on Magnitude (and Rupture areas extent). During the selection of rupture area barycenter, with this choice, the nodes closer than 0.25 times the Width to the mesh edge are discarded.
	"minimum_interdistance": 0.1,           # With this choices, the selected rupture barycenters are distant from each other more than 0.1 times the Length. This will avoid to have very similar rupture areas and reduce the number of scenarios at largest magnitude bins (see Scala et al. 2020) 
	}
	}

 ## 2.2 scaling_relationship.json

The magnitude bins and the rupture geometries (according to the selected scaling laws) are set in the input file [scaling_relationship.json](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/config_files/Parameters/scaling_relationship.json) contained in the [config_files/Parameters](https://github.com/antonioscalaunina/pyANTI-FASc/tree/main/config_files/Parameters) folder. The file actually used for the run must be always named scaling_relationship.json. The default available file is set to run this example where we use a selection similar to that one proposed in the framework of the project TSUMAPS-NEAM (see Basili et al. 2021) using the Strasser et al. (2010) and the Murotani et al.(2013) scaling relationships. However, it can be easily modified to run the software with different magnitude binnings and scaling laws. In the example below the structure of this file, look carefully at the comments beside to properly set the values:

    { 
    "Magnitude_bins": {                                     # Within this section the number of magnitude bins and the magnitude bins are defined
    "number_bins" : 32, 
    "Magnitude": [6.0000, 6.5000, 6.8012, 7.0737, 7.3203, 7.5435, 7.7453, 7.9280, 8.0933,
              8.2429, 8.3782, 8.5007, 8.6115, 8.7118, 8.8025, 8.8846, 8.9588, 9.0260, 
	      9.0869, 9.1419, 9.1917, 9.2367, 9.2775, 9.3144, 9.3478, 9.3780, 9.4053, 
	      9.4300, 9.4524 , 9.4727, 9.4910,9.5075]
    },

    "Scaling_law": { "number": 2,                            # Here we declare the number of different scaling laws used in the code        
    "name" : ["Murotani", "Strasser"],                       # Names of scaling laws (must be consistent with the selected number
    "Area": [156.233, 494.051, 988.488, 1851.277, 3266.416, 5460.991, 8691.034, 13236.448,                         #Values of the area. They must be "number_bins" * "number" (of Scaling law). In this case the first 32 refer to Murotani scaling law, while the remaining one to the Strasser scaling
	          19367.407, 27332.004, 37322.564, 49484.760, 63866.258,  80458.414, 99145.149,
	          119776.459, 142092.603, 165871.385, 190821.652, 216599.421, 242912.949, 
	          269467.272, 295988.071, 322227.910, 347970.147, 373030.811, 397258.801, 
	          420534.777, 442769.108, 463899.186, 483886.383, 502712.864,
        172.187, 515.229, 997.108, 1812.018, 3111.183, 5074.719, 7898.154, 11788.431, 
	          16936.419, 23509.360, 31626.155, 41368.179, 52740.956, 65710.323, 80164.115, 
	      	  95970.819, 112921.750, 130843.454, 149515.800, 168684.550, 188138.900, 207668.823, 
		  227081.948, 246207.279, 264897.303, 283028.731, 300502.136, 317240.766, 333188.779, 
		  348309.117, 362581.201, 375998.580],
    "Length": [10.2770, 20.7196, 31.6095, 46.3198, 65.4563, 89.5107, 118.7881, 153.4749, 193.5108,                #Values of the length. They must be "number_bins" * "number" (of Scaling law). In this case the first 32 refer to Murotani scaling law, while the remaining one to the Strasser scaling. The Width W will be computed for each bin as Area/Length
                    238.6781, 288.5428, 342.6192, 400.2115, 460.6512, 523.1286, 586.9580,
                    651.3207, 715.6812, 779.437, 841.965, 902.856, 961.739, 
		    1018.321, 1072.384, 1123.771, 1172.388, 1218.190, 
		    1261.171, 1301.370, 1338.847, 1373.687, 1405.992,
           10.789, 21.159, 31.747, 45.826, 63.882, 86.287, 113.240, 144.837, 180.959, 
                    221.359, 265.612, 313.263, 363.687, 416.297, 470.395, 525.401, 580.628, 
		    635.638,689.940, 743.025, 794.570, 844.286, 891.947, 937.387, 980.496, 
		    1021.208, 1059.502, 1095.388, 1128.905, 1160.116, 1189.101, 1215.950]
 
    }
    }

# 3 Run pyANTI-FASc

Once the mesh is selected and the other configuration paramters are set through the input files, the whole process can be launched simply running the following command:

	python antifasc_main.py

 The output on the screen will allow the user to follow the different steps of the running and that everything is working. Below some details.

 The software reads the input files and let the user know that the selected mesh discretization (nodes and cells) has been found. The mesh is then written in an inp format. 

 	reading input.json file
	reading scaling_relationship.json file
	Great! You already have nodes and cells, I'm just writing

 Subsequently, the software starts the barycenter selection according to the selected application. As you can see in the example, only magnitude bins around the real magnitude of the event is used accordingly to the selected application (PTF). The output on the screen confirms that the two selected scaling laws have been used

	Barycenter selection
	Magnitude bin # 15 - Mw=8.8846
	Magnitude bin # 16 - Mw=8.9588
	Magnitude bin # 17 - Mw=9.0260
	Magnitude bin # 18 - Mw=9.0869
	Magnitude bin # 19 - Mw=9.1419
	Barycenter selection (PTF)
	Magnitude bin # 15 - Mw=8.8846
	Magnitude bin # 16 - Mw=8.9588
	Magnitude bin # 17 - Mw=9.0260
	Magnitude bin # 18 - Mw=9.0869
	Magnitude bin # 19 - Mw=9.1419
	Mw: [8.8846 8.9588 9.026  9.0869 9.1419]
	Scaling names: ['Murotani', 'Strasser']

After that, the rupture areas computation starts, for each bin of magnitude and each scaling law. For each of these classes the output on the screen indicates how many rupture areas have been computed. Finally the rupture areas are written in temporary output files that will be then  used as input for the slip distribution defintions
	
 	Rupturing area computation
	Magnitude bin # 15 - Mw=8.8846
	Magnitude bin # 16 - Mw=8.9588
	Magnitude bin # 17 - Mw=9.0260
	Magnitude bin # 18 - Mw=9.0869
	Magnitude bin # 19 - Mw=9.1419
	Rupturing areas computed!
	Mw=8.8846, Name scaling: Murotani, N=36, N_all=36
	Mw=8.8846, Name scaling: Strasser, N=40, N_all=40
	Mw=8.9588, Name scaling: Murotani, N=31, N_all=31
	Mw=8.9588, Name scaling: Strasser, N=36, N_all=36
	Mw=9.026, Name scaling: Murotani, N=26, N_all=26
	Mw=9.026, Name scaling: Strasser, N=33, N_all=33
	Mw=9.0869, Name scaling: Murotani, N=23, N_all=23
	Mw=9.0869, Name scaling: Strasser, N=29, N_all=29
	Mw=9.1419, Name scaling: Murotani, N=16, N_all=16
	Mw=9.1419, Name scaling: Strasser, N=25, N_all=25
	Writing Output
	Magnitude bin # 15 - Mw=8.8846
	Magnitude bin # 16 - Mw=8.9588
	Magnitude bin # 17 - Mw=9.0260
	Magnitude bin # 18 - Mw=9.0869
	Magnitude bin # 19 - Mw=9.1419

 The slip distributions will be then computed and the screen output lets the user know within each class the software is working

 	Computing slip distributions for the homogeneous and variable rigidity cases
	
 	/home/scala/pyANTI-FASc/Tohoku_test_M90_E14237_N3832_slip_KuJ/homogeneous_mu/8_8846/Murotani
 	starting ...
 	Number of scenarios is         180
	/home/scala/pyANTI-FASc/Tohoku_test_M90_E14237_N3832_slip_KuJ/homogeneous_mu/8_8846/Strasser
 	starting ...
 	Number of scenarios is         200
	/home/scala/pyANTI-FASc/Tohoku_test_M90_E14237_N3832_slip_KuJ/homogeneous_mu/8_9588/Murotani
 	starting ...
 	Number of scenarios is         155
	/home/scala/pyANTI-FASc/Tohoku_test_M90_E14237_N3832_slip_KuJ/homogeneous_mu/8_9588/Strasser
 	starting ...
 	Number of scenarios is         180
	/home/scala/pyANTI-FASc/Tohoku_test_M90_E14237_N3832_slip_KuJ/homogeneous_mu/9_0260/Murotani
 	starting ...
 	Number of scenarios is         130
	/home/scala/pyANTI-FASc/Tohoku_test_M90_E14237_N3832_slip_KuJ/homogeneous_mu/9_0260/Strasser
 	starting ...
 	Number of scenarios is         165
	.......

 It is worth to highlight that the number of scenarios, that is the number of slip distributions is always 5 times the number of selected areas as set in the input configuration file.

 When the process is over, The output will be finally organized as shown in the following tree:

 	output
	├── Tohoku_test_M90_E14237_N3832_slip_KuJ
	│   ├── homogeneous_mu
	│   │   ├── 8_8846
	│   │   │   ├── Murotani
	│   │   │   │   ├── Slip4HySea00004_001.dat
	│   │   │   │   ├── Slip4HySea00004_002.dat
	│   │   │   │   ├── Slip4HySea00004_003.dat
	│   │   │   │   ├── Slip4HySea00004_004.dat
	│   │   │   │   ├── Slip4HySea00004_005.dat
	│   │   │   │   ├── Slip4HySea00007_001.dat
	│   │   │   │   ├── Slip4HySea00007_002.dat
	│   │   │   │   ├── Slip4HySea00007_003.dat
	│   │   │   │   ├── Slip4HySea00007_004.dat
	│   │   │   │   ├── Slip4HySea00007_005.dat
	│   │   │   │   ├── Slip4HySea00008_001.dat
	│   │   │   │   ├── Slip4HySea00008_002.dat
	│   │   │   │   ├── Slip4HySea00008_003.dat
	│   │   │   │   ├── Slip4HySea00008_004.dat
	│   │   │   │   ├── Slip4HySea00008_005.dat
	│   │   │   │   ├── Slip4HySea00023_001.dat
	│   │   │   │   ├── Slip4HySea00023_002.dat
	│   │   │   │   ├── Slip4HySea00023_003.dat
	│   │   │   │   ├── Slip4HySea00023_004.dat
	│   │   │   │   ├── Slip4HySea00023_005.dat
	│   │   │   │   ├── Slip4HySea00027_001.dat
	│   │   │   │   ├── Slip4HySea00027_002.dat
	│   │   │   │   ├── Slip4HySea00027_003.dat
	│   │   │   │   ├── Slip4HySea00027_004.dat
	│   │   │   │   ├── Slip4HySea00027_005.dat
 	.............................................

The Slip4HySea* files are in the standard format used as input by the software [Tsunami-HySea](https://edanya.uma.es/hysea/):

	 LON1     LAT1    DEPTH1(km)      LON2    LAT2    DEPTH2(km)      LON3    LAT3    DEPTH3(km)      RAKE    SLIP(m)
	  142.209106   36.309288    7.494581  142.289078   36.401661    7.454205  142.143066   36.442127   10.402809   90.000000    7.256678
	  142.209106   36.309288    7.494581  142.336136   36.256641    5.620936  142.289078   36.401661    7.454205   90.000000    7.095560
	  142.289078   36.401661    7.454205  142.223038   36.534504   10.362430  142.143066   36.442127   10.402809   90.000000    6.840256
	  142.209106   36.309288    7.494581  142.143066   36.442127   10.402809  142.074631   36.344494   10.356820   90.000000    7.429152
	  142.209106   36.309288    7.494581  142.275635   36.170479    5.620936  142.336136   36.256641    5.620936   90.000000    6.319403
	  142.336136   36.256641    5.620936  142.416107   36.349018    5.580560  142.289078   36.401661    7.454205   90.000000    6.555849
	  142.369064   36.494038    7.413829  142.223038   36.534504   10.362430  142.289078   36.401661    7.454205   90.000000    6.771578
	  142.223038   36.534504   10.362430  142.077026   36.574970   13.311030  142.143066   36.442127   10.402809   90.000000    5.898129
	  142.074631   36.344494   10.356820  142.143066   36.442127   10.402809  142.008591   36.477337   13.265040   90.000000    6.484488
	  142.209106   36.309288    7.494581  142.074631   36.344494   10.356820  142.143524   36.212208    7.409739   90.000000    7.508512
	  142.209106   36.309288    7.494581  142.143524   36.212208    7.409739  142.275635   36.170479    5.620936   90.000000    6.776600
	  142.275635   36.170479    5.620936  142.402664   36.117832    3.747290  142.336136   36.256641    5.620936   90.000000    5.396324
	  142.463165   36.203999    3.747290  142.416107   36.349018    5.580560  142.336136   36.256641    5.620936   90.000000    6.399800
	  142.416107   36.349018    5.580560  142.369064   36.494038    7.413829  142.289078   36.401661    7.454205   90.000000    6.580969
	  142.369064   36.494038    7.413829  142.303024   36.626881   10.322050  142.223038   36.534504   10.362430   90.000000    6.556834
   .........................................................................................................................................

## 4 Post-process

 The slip distribution can be easily plotted by simple personal scripts. In the folder [utils](https://github.com/antonioscalaunina/pyANTI-FASc/tree/main/utils) there is the script [plot_slip_distribution.py](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/utils/plot_slip_distribution.py). It is launched with the following command:
 
	python plot_slip_distribution.py

 This script will ask which class the user wants to plot as shown in the example below:

 	Current folder is '../output/'

	There is only one event directory

	../output/Tohoku_test_M90_E14237_N3832_slip_KuJ/
	Current folder is '../output/Tohoku_test_M90_E14237_N3832_slip_KuJ/'

	Choose your rigidity distribution directory between:
	1. homogeneous_mu/
	2. variable_mu/
	Insert a number between 1 and 2:

	2
	Current folder is '../output/Tohoku_test_M90_E14237_N3832_slip_KuJ/variable_mu/'

	Choose your magnitude directory between:
	1. 8_8846/
	2. 8_9588/
	3. 9_0260/
	4. 9_0869/
	5. 9_1419/
	Insert a number between 1 and 5:
	
	3
	Current folder is '../output/Tohoku_test_M90_E14237_N3832_slip_KuJ/variable_mu/9_0260/'
	
	Choose your scaling law directory between:
	1. Murotani/
	2. Strasser/
	Insert a number between 1 and 2:
	
	1
	Number of files:  130
	[##############################]  100%

 In this example we have finally selected the folder */output/Tohoku_test_M90_E14237_N3832_slip_KuJ/variable_mu/9_0260/*. Within that folder, for each Slip4Hysea*.dat file, two new files will be produced: one is the slip distribution in a standard geoJSON format (you might find an example [here](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/utils/Slip4HySea00007_004.json)) and the other an interactive maps in HTML format (example [here](https://antonioscalaunina.github.io/pyANTI-FASc/utils/Slip4HySea00007_004.html)).

 
The test-case shown in this example is also run through the Jupyter Notebook [antifasc_main.ipynb](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/bin/antifasc_main.ipynb). Within the Jupyter Notebook some intermediate plots are shown to better outline some of the steps of the process, such as the barycenter selection, and the rupture area computation.


   
