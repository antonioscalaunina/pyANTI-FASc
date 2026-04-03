# EXAMPLE 1 Simulation in the Mediterranean (still incomplete)
## from EFSM20 services

In this brief guide, a practical example to run a pyANTI-FASc application is shown. 
In this example, a set of slip distributions simulating an earthquake offshore Sicily, southern Italy, is generated.
The test-case shown in this example is also run through the Jupyter Notebook [antifasc_main.ipynb](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/bin/antifasc_main.ipynb). Within the Jupyter Notebook some intermediate plots are shown to better describe some of the steps of the process like the barycenter selection and the rupture area computation.


# 1 - Mesh

pyANTI-FASc is able to generate outputs using meshes made available by the European Database of Seismogenic Faults services ([EFSM20](https://seismofaults.eu/component/tags/tag/efsm20)).

In particular, it supports the meshes provided by the *EFSM20 Meshes* service, whose details are available [here](https://seismofaults.eu/services/efsm20-services). At this [link](https://services.seismofaults.eu/EFSM20_Meshes/ows?service=WFS&request=GetCapabilities), you can access the XML file required to load the database into QGIS software.

If you use this service, please cite: [https://doi.org/10.13127/efsm20/meshes](https://doi.org/10.13127/efsm20/meshes).

Figure 1 shows a screenshot of some of the meshed faults available within the EFSM20 database.


![Screenshot of EFSM mesh database imported within QGis](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/utils/sz_slabs/Image_qgis.jpg)
*Figure 1 - Screenshot of EFSM mesh database imported within QGis*


# 2 - Input files

## 2.1 input.json

This file [input_Sicily.json](https://github.com/antonioscalaunina/pyANTI-FASc/tree/main/config_files/Parameters/input_Sicily.json) contained in the [config_files/Parameters](https://github.com/antonioscalaunina/pyANTI-FASc/tree/main/config_files/Parameters) folder manages the main configuration parameters to run this examples. The file actually used for the run must be always named **input.json**. The default available file is set to run the Tohoku example, but it is sufficient to copy this `input_Sicily.json` into `input.json` to run this example

Here below the important settings to be managed by the user are shown. **Look carefully to the comments besides and below the parameters**. The parameters not shown in this example might be left unmodified. Their use and functionality will be fixed in next releases and better described in [Wiki Documentation](https://github.com/antonioscalaunina/pyANTI-FASc/wiki) currently under construction.

	{"zone_name": "ITCF00G",      # Name of the precomputed mesh to be used. The mesh downloaded from the ESFM20 service must be saved in geojson file as within the `utils/sz_slabs/' folder as shown [here](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/utils/sz_slabs/ITCF00G_mesh.json)
	"Merc_zone": 33,                  # Mercator zone for the selected slab. See the slab database and use the correct Mercator zone 
    "acronym": "ITC",		  # 3 digit acronym that is used for that slab. It can be arbitrarily chosen by the user (but must have 3 digits!). You might find suggestions into the slab database.
	
	"mesh_gen": 1,           
	"rake": -90,

	# If the "mesh_gen" option is set to 1, a GeoJSON file containing the mesh is expected.
	# The file name must be consistent with the "zone_name" (e.g., zone_name_mesh.json),
	# and it should be located in utils/sz_slabs 
	# (see example:https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/utils/sz_slabs/ITCF00G_mesh.json in the repository.
	# If the GeoJSON file is not found or the name is incorrect, the run will stop and
	# the following error message will be shown to the user:
	# ERROR: Mesh in GeoJSON format does not exist! Please check the options in input.json and the zone/file names.

	# The "rake" option allows the user to specify a rake value.
	# If not provided, a default value of rake = 90° is assigned to the entire mesh.
	# Alternatively, this field can contain the path to a CSV file where rake values are defined.
	# In this case, the number of rake values must match the number of mesh cells (and defined in the same order).
	# If this condition is not met, the default rake = 90° will be assigned to the entire mesh.
	
We propose a run in "Hazard" mode (see below). All the possible slip distributions (with an optimized number of rupture areas decreasing with magnitude) is computed. The section event can be hence left unmodified
       	
	"Configure": {
	"application": "Hazard",                     # This application restricts the computed scenarios to a range of magnitude and location around predefined values
	"shape": "Rectangle",                     # This choice allows to compute scenarios with aspect ratio L/W preserved as prescribed by the selected scaling law. The other possible choice is "Circle". More details soon in the Wiki Documentation
	"numb_stoch": 5,                          # Number of stochastic slip for each rupture areas
	"variable_mu": 1,                         # 1 means that also the distributions with variable rigidity will be computed. 0 for computing only the case with homogeneous rigidity
	"coupling_shallow_limit":1.0,             # Shallow limit for the area where the seismic coupling is expected to decrease
	"coupling_deep_limit":55.0,		  # Deep limit for the area where the seismic coupling is expected to decrease


 
	
	"minimum_bnd_distance": 0.1,           # This option (as well as the next one) is used to limit the number of rupture areas dependending on Magnitude (and Rupture areas extent). During the selection of rupture area barycenter, with this choice, the nodes closer than 0.25 times the Width to the mesh edge are discarded.
	"minimum_interdistance": 0.2,           # With this choices, the distance between the selected rupture barycenters will be more than 0.1 times the Length. This will avoid to have very similar rupture areas and reduce the number of scenarios at largest magnitude bins (see Scala et al. 2020) 
	}
	}


 ## 2.2 scaling_relationship.json

 The magnitude bins and the rupture geometries (defined by the selected scaling laws) for this application are set in the input file [scaling_relationship_WC.json](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/config_files/Parameters/scaling_relationship_WC.json) contained in the [config_files/Parameters](https://github.com/antonioscalaunina/pyANTI-FASc/tree/main/config_files/Parameters) folder. In this file, the [Wells & Coppersmith 1994](https://doi.org/10.1785/BSSA0840040974) scaling relationship for normal faulting is implemented.
The file actually used for the run must be always named **scaling_relationship.json**. The default available file is set to run the [Tohoku example](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/Example1_Tohoku.md), but it is sufficient to copy this `scaling_relationship_WC.json` into `scaling_relathionship.json` to run this example**look carefully at the comments beside and below to properly set the values**:

    { 
    "Magnitude_bins": {                                     # Within this section the number of magnitude bins and the magnitude bins are defined
    "number_bins" : 21, 
	"Magnitude": [5.5,5.6,5.7,5.8,5.9,6.0,6.1,6.2,6.3,6.4,6.5,6.6,6.7,6.8,6.9,7.0,7.1,7.2,7.3,7.4,7.5]
	},

	"Scaling_law": { "number": 1,      # Here we declare the number of different scaling laws used in the code 
	"name" : ["WellsCopp1994_Normal"],        # Names of scaling laws (must be consistent with the set number in the parameter "Scaling_law"
	
	"Area":  [43.6516,   52.7230,   63.6796,   76.9130,   92.8966,  112.2018,  135.5189,  163.6817,  
	197.6970,  238.7811,  288.4032,  348.3373,  420.7266,  508.1594,  	613.7620,  741.3102,  
	895.3648, 1081.4340, 1306.1709, 1577.6113, 1905.4607],  
	"Length": [7.4131,    8.3176,    9.3325,   10.4713,   11.7490,   13.1826,   14.7911,   16.5959,  
	18.6209,   20.8930,   23.4423,   26.3027,   29.5121,   33.1131,   37.1535,   41.6869,   46.7735,  
	52.4807,   58.8844,   66.0693,   74.1310] 
	}
	}
	# Number of area and length values must be "number_bins" * "number" (of Scaling law). 

# 3 Run pyANTI-FASc

Once the mesh is selected and the other configuration parameters are set through the described input files, the whole process can be launched simply running the following commands:

	conda activate antifasc
	cd bin
	python antifasc_main.py

 The output on the screen will allow the user to follow the different steps of the running and that everything is working. Below some details.

 The software reads the input files and let the user know that the selected mesh discretization (nodes and cells) has been found and will be used for the run
 	
  	reading input.json file
	reading scaling_relationship.json file
	../utils/sz_slabs/ITCF00G_mesh.json
	Great! I really love to create mesh files from GeoJSON format

	Files ../utils/sz_slabs/ITCF00G/subfaults/ITCF00G_mesh_nodes.dat and ../utils/sz_slabs/ITCF00G/subfaults/ITCF00G_mesh_faces.dat created successfully inside subfaults!
This last message inform the author that from now on, the same mesh can be used setting `"mesh_gen":0,` in the `input.json` file, since the new mesh is now part of the `utils/sz_slabs` databese 

 Subsequently, the software performs the barycenter selection according to the selected application. As you can see in the example, all the magnitude bins in `scaling_relationship.json` file are selected according to the application `Hazard`. The output on the screen confirms that the selected scaling law has been used

	Barycenter selection
	Magnitude bin # 0 - Mw=5.5000
	Magnitude bin # 1 - Mw=5.6000
	Magnitude bin # 2 - Mw=5.7000
	Magnitude bin # 3 - Mw=5.8000
	Magnitude bin # 4 - Mw=5.9000
	Magnitude bin # 5 - Mw=6.0000
	Magnitude bin # 6 - Mw=6.1000
	Magnitude bin # 7 - Mw=6.2000
	Magnitude bin # 8 - Mw=6.3000
	Magnitude bin # 9 - Mw=6.4000
	Magnitude bin # 10 - Mw=6.5000
	Magnitude bin # 11 - Mw=6.6000
	Magnitude bin # 12 - Mw=6.7000
	Magnitude bin # 13 - Mw=6.8000
	Magnitude bin # 14 - Mw=6.9000
	Magnitude bin # 15 - Mw=7.0000
	Magnitude bin # 16 - Mw=7.1000
	Magnitude bin # 17 - Mw=7.2000
	Magnitude bin # 18 - Mw=7.3000
	Magnitude bin # 19 - Mw=7.4000
	Magnitude bin # 20 - Mw=7.5000
	Mw: [5.5 5.6 5.7 5.8 5.9 6.  6.1 6.2 6.3 6.4 6.5 6.6 6.7 6.8 6.9 7.  7.1 7.2
	 7.3 7.4 7.5]
	Scaling names: ['WellsCopp1994_Normal']

After that, the rupture areas computation is performed (a waiting bar informs the user about the advancement), for each bin of magnitude and for the selected scaling law. For each of these classes, the output on the screen indicates how many rupture areas have been computed. Finally the rupture areas are written in temporary output files that will be then  used as input for the slip distribution computation
	
	Computing Rupturing areas: 100%|████████████████████████████████████████████████████████| 21/21 [00:01<00:00, 11.74it/s]
	Mw=5.5, Name scaling: WellsCopp1994_Normal, N=743, N_all=743
	Mw=5.6, Name scaling: WellsCopp1994_Normal, N=743, N_all=743
	Mw=5.7, Name scaling: WellsCopp1994_Normal, N=743, N_all=743
	Mw=5.8, Name scaling: WellsCopp1994_Normal, N=743, N_all=743
	Mw=5.9, Name scaling: WellsCopp1994_Normal, N=743, N_all=743
	Mw=6.0, Name scaling: WellsCopp1994_Normal, N=737, N_all=737
	Mw=6.1, Name scaling: WellsCopp1994_Normal, N=737, N_all=737
	Mw=6.2, Name scaling: WellsCopp1994_Normal, N=735, N_all=735
	Mw=6.3, Name scaling: WellsCopp1994_Normal, N=735, N_all=735
	Mw=6.4, Name scaling: WellsCopp1994_Normal, N=312, N_all=312
	Mw=6.5, Name scaling: WellsCopp1994_Normal, N=312, N_all=312
	Mw=6.6, Name scaling: WellsCopp1994_Normal, N=199, N_all=199
	Mw=6.7, Name scaling: WellsCopp1994_Normal, N=149, N_all=149
	Mw=6.8, Name scaling: WellsCopp1994_Normal, N=128, N_all=128
	Mw=6.9, Name scaling: WellsCopp1994_Normal, N=100, N_all=100
	Mw=7.0, Name scaling: WellsCopp1994_Normal, N=82, N_all=82
	Mw=7.1, Name scaling: WellsCopp1994_Normal, N=69, N_all=69
	Mw=7.2, Name scaling: WellsCopp1994_Normal, N=54, N_all=54
	Mw=7.3, Name scaling: WellsCopp1994_Normal, N=44, N_all=44
	Mw=7.4, Name scaling: WellsCopp1994_Normal, N=37, N_all=37
	Mw=7.5, Name scaling: WellsCopp1994_Normal, N=32, N_all=32
	Writing Output
	Magnitude bin # 0 - Mw=5.5000
	Magnitude bin # 1 - Mw=5.6000
	Magnitude bin # 2 - Mw=5.7000
	Magnitude bin # 3 - Mw=5.8000
	Magnitude bin # 4 - Mw=5.9000
	Magnitude bin # 5 - Mw=6.0000
	Magnitude bin # 6 - Mw=6.1000
	Magnitude bin # 7 - Mw=6.2000
	Magnitude bin # 8 - Mw=6.3000
	Magnitude bin # 9 - Mw=6.4000
	Magnitude bin # 10 - Mw=6.5000
	Magnitude bin # 11 - Mw=6.6000
	Magnitude bin # 12 - Mw=6.7000
	Magnitude bin # 13 - Mw=6.8000
	Magnitude bin # 14 - Mw=6.9000
	Magnitude bin # 15 - Mw=7.0000
	Magnitude bin # 16 - Mw=7.1000
	Magnitude bin # 17 - Mw=7.2000
	Magnitude bin # 18 - Mw=7.3000
	Magnitude bin # 19 - Mw=7.4000
	Magnitude bin # 20 - Mw=7.5000

 The slip distributions is then finally computed and the screen standard output lets the user know within each class the software is working

 	Computing slip distributions for the homogeneous and variable rigidity cases
	
 	/mnt/c/Users/ascal/Downloads/pyANTI-FASc/Tohoku_test_M90_E14237_N3832_slip_KuJ/homogeneous_mu/8_9588/Murotani
	 starting ...
	 Number of scenarios is          50
	/mnt/c/Users/ascal/Downloads/pyANTI-FASc/Tohoku_test_M90_E14237_N3832_slip_KuJ/homogeneous_mu/8_9588/Strasser
	 starting ...
	 Number of scenarios is          55
	/mnt/c/Users/ascal/Downloads/pyANTI-FASc/Tohoku_test_M90_E14237_N3832_slip_KuJ/homogeneous_mu/9_0260/Murotani
	 starting ...
	 Number of scenarios is          35
	/mnt/c/Users/ascal/Downloads/pyANTI-FASc/Tohoku_test_M90_E14237_N3832_slip_KuJ/homogeneous_mu/9_0260/Strasser
	 starting ...
	 Number of scenarios is          50
	/mnt/c/Users/ascal/Downloads/pyANTI-FASc/Tohoku_test_M90_E14237_N3832_slip_KuJ/homogeneous_mu/9_0869/Murotani
	 starting ...
	 Number of scenarios is          35
	/mnt/c/Users/ascal/Downloads/pyANTI-FASc/Tohoku_test_M90_E14237_N3832_slip_KuJ/homogeneous_mu/9_0869/Strasser
	 starting ...
	 Number of scenarios is          35
	/mnt/c/Users/ascal/Downloads/pyANTI-FASc/Tohoku_test_M90_E14237_N3832_slip_KuJ/variable_mu/8_9588/Murotani
	 starting ...
	 Number of scenarios is          50
	/mnt/c/Users/ascal/Downloads/pyANTI-FASc/Tohoku_test_M90_E14237_N3832_slip_KuJ/variable_mu/8_9588/Strasser
	 starting ...
	 Number of scenarios is          55
	/mnt/c/Users/ascal/Downloads/pyANTI-FASc/Tohoku_test_M90_E14237_N3832_slip_KuJ/variable_mu/9_0260/Murotani
	 starting ...
	 Number of scenarios is          35
	/mnt/c/Users/ascal/Downloads/pyANTI-FASc/Tohoku_test_M90_E14237_N3832_slip_KuJ/variable_mu/9_0260/Strasser
	 starting ...
	 Number of scenarios is          50
	/mnt/c/Users/ascal/Downloads/pyANTI-FASc/Tohoku_test_M90_E14237_N3832_slip_KuJ/variable_mu/9_0869/Murotani
	 starting ...
	 Number of scenarios is          35
	/mnt/c/Users/ascal/Downloads/pyANTI-FASc/Tohoku_test_M90_E14237_N3832_slip_KuJ/variable_mu/9_0869/Strasser
	 starting ...
	 Number of scenarios is          35

 It is worth to highlight that the number of scenarios, that is the number of slip distributions is always 5 times the number of selected areas as set in the **input.json** configuration file.

 When the process is over, The output will be finally organized as shown in the following tree:

 	output
	├── Tohoku_test_M90_E14237_N3832_slip_KuJ
	│   ├── homogeneous_mu
	│   │   ├── 8_9588
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

The Slip4HySea* files are in the standard format used as input by the software [Tsunami-HySea](https://edanya.uma.es/hysea/) which is one of the most widely used tsunami simulators within the community.

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

 The slip distributions can be easily plotted by simple personal scripts. In the folder [utils](https://github.com/antonioscalaunina/pyANTI-FASc/tree/main/utils) there is the script [plot_slip_distribution.py](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/utils/plot_slip_distribution.py). It might be run (still in the *antifasc* Conda environment) with the following command:
 
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
	1. 8_9588/
	2. 9_0260/
	3. 9_0869/
	Insert a number between 1 and 3:
	
	2
	Current folder is '../output/Tohoku_test_M90_E14237_N3832_slip_KuJ/variable_mu/9_0260/'
	
	Choose your scaling law directory between:
	1. Murotani/
	2. Strasser/
	Insert a number between 1 and 2:
	
	2
	Number of files:  50
	[##############################]  100%

 In this example we have finally selected the folder */output/Tohoku_test_M90_E14237_N3832_slip_KuJ/variable_mu/9_0260/*. Within that folder, for each Slip4Hysea*.dat file, two new files will be produced: one is the slip distribution in a standard geoJSON format (you might find an example [here](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/utils/Slip4HySea00004_002.json)) while the other one is an interactive maps in HTML format (example [here](https://antonioscalaunina.github.io/pyANTI-FASc/utils/Slip4HySea00004_002.html)).

As already outlined in the [README](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/README.md), this example can be run through the Jupyter Notebook available [here](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/bin/antifasc_main.ipynb). For example, you can install an IDE like [Visual Studio Code](https://code.visualstudio.com/download)  to run the Jupyter Notebook on Ubuntu, macOS, or Windows. 
**IMPORTANT**: To ensure everything works correctly, it's essential to properly install all the needed components as outlined in the README.


