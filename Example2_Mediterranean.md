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


![Screenshot of EFSM mesh database imported to QGis](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/utils/sz_slabs/Image_qgis.jpg)
*Figure 1 - Screenshot of EFSM mesh database imported to QGis*


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
The file actually used for the run must be always named **scaling_relationship.json**. The default available file is set to run the [Tohoku example](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/Example1_Tohoku.md), but it is sufficient to copy this `scaling_relationship_WC.json` into `scaling_relathionship.json` to run this example **look carefully at the comments beside and below to properly set the values**:

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

 The slip distributions is then finally computed and the screen standard output lets the user know within each class the software is working. The last line informs the user about the overall running time

 	Computing slip distributions for the homogeneous and variable rigidity cases
	/home/scala/pyANTI-FASc/ITCF00G_Hazard_slip_ITC/homogeneous_mu/5_5000/WellsCopp1994_Normal
	 starting ...
	 Number of scenarios is        3715
	 scenario #          250
	 scenario #          500
	 scenario #          750
	 scenario #         1000
	 scenario #         1250
	 scenario #         1500
	 scenario #         1750
	 scenario #         2000
	 scenario #         2250
	 scenario #         2500
	 scenario #         2750
	 scenario #         3000
	 scenario #         3250
	 scenario #         3500
	/home/scala/pyANTI-FASc/ITCF00G_Hazard_slip_ITC/homogeneous_mu/5_6000/WellsCopp1994_Normal
	 starting ...
	 Number of scenarios is        3715
	 scenario #          250
	 scenario #          500
	 scenario #          750
	 scenario #         1000
	 scenario #         1250
	 scenario #         1500
	 scenario #         1750
	 scenario #         2000
	 scenario #         2250
	 scenario #         2500
	 scenario #         2750
	 scenario #         3000
	 scenario #         3250
	 scenario #         3500
	 /home/scala/pyANTI-FASc/ITCF00G_Hazard_slip_ITC/homogeneous_mu/5_7000/WellsCopp1994_Normal
	..........................................................................................ù
	 scenario #          250
	 scenario #          500
	 scenario #          750
	 scenario #         1000
	 scenario #         1250
	 scenario #         1500
	/home/scala/pyANTI-FASc/ITCF00G_Hazard_slip_ITC/variable_mu/6_6000/WellsCopp1994_Normal
	 starting ...
	 Number of scenarios is         995
	 scenario #          250
	 scenario #          500
	 scenario #          750
	/home/scala/pyANTI-FASc/ITCF00G_Hazard_slip_ITC/variable_mu/6_7000/WellsCopp1994_Normal
	 starting ...
	 Number of scenarios is         745
	 scenario #          250
	 scenario #          500
	/home/scala/pyANTI-FASc/ITCF00G_Hazard_slip_ITC/variable_mu/6_8000/WellsCopp1994_Normal
	 starting ...
	 Number of scenarios is         640
	 scenario #          250
	 scenario #          500
	/home/scala/pyANTI-FASc/ITCF00G_Hazard_slip_ITC/variable_mu/6_9000/WellsCopp1994_Normal
	 starting ...
	 Number of scenarios is         500
	 scenario #          250
	 scenario #          500
	/home/scala/pyANTI-FASc/ITCF00G_Hazard_slip_ITC/variable_mu/7_0000/WellsCopp1994_Normal
	 starting ...
	 Number of scenarios is         410
	 scenario #          250
	/home/scala/pyANTI-FASc/ITCF00G_Hazard_slip_ITC/variable_mu/7_1000/WellsCopp1994_Normal
	 starting ...
	 Number of scenarios is         345
	 scenario #          250
	/home/scala/pyANTI-FASc/ITCF00G_Hazard_slip_ITC/variable_mu/7_2000/WellsCopp1994_Normal
	 starting ...
	 Number of scenarios is         270
	 scenario #          250
	/home/scala/pyANTI-FASc/ITCF00G_Hazard_slip_ITC/variable_mu/7_3000/WellsCopp1994_Normal
	 starting ...
	 Number of scenarios is         220
	/home/scala/pyANTI-FASc/ITCF00G_Hazard_slip_ITC/variable_mu/7_4000/WellsCopp1994_Normal
	 starting ...
	 Number of scenarios is         185
	/home/scala/pyANTI-FASc/ITCF00G_Hazard_slip_ITC/variable_mu/7_5000/WellsCopp1994_Normal
	 starting ...
	 Number of scenarios is         160
	316.18324065208435
 It is worth to highlight that the number of scenarios, that is the number of slip distributions is always 5 times the number of selected areas as set in the `input.json` configuration file.

 When the process is over, The output will be finally organized as shown in the following tree:

 	output/ITCF00G_Hazard_slip_ITC/
	├── homogeneous_mu
	│   ├── 5_5000
	│   │   └── WellsCopp1994_Normal
	│   │       ├── output_file.log
	│   │       ├── Slip4HySea00001_000.dat
	│   │       ├── Slip4HySea00002_000.dat
	│   │       ├── Slip4HySea00003_000.dat
	│   │       ├── Slip4HySea00004_000.dat
	│   │       ├── Slip4HySea00005_000.dat
	│   │       ├── Slip4HySea00006_000.dat
	│   │       ├── Slip4HySea00007_000.dat
	│   │       ├── Slip4HySea00008_000.dat
	│   │       ├── Slip4HySea00009_000.dat
	│   │       ├── Slip4HySea00010_000.dat
	│   │       ├── Slip4HySea00011_000.dat
	│   │       ├── Slip4HySea00012_000.dat
	│   │       ├── Slip4HySea00013_000.dat
	│   │       ├── Slip4HySea00014_000.dat
	│   │       ├── Slip4HySea00015_000.dat
	│   │       ├── Slip4HySea00016_000.dat
	│   │       ├── Slip4HySea00017_000.dat
	│   │       ├── Slip4HySea00018_000.dat
	│   │       ├── Slip4HySea00019_000.dat
	│   │       ├── Slip4HySea00020_000.dat
	│   │       ├── Slip4HySea00021_000.dat
	│   │       ├── Slip4HySea00022_000.dat
	│   │       ├── Slip4HySea00023_000.dat
	│   │       ├── Slip4HySea00024_000.dat
	│   │       ├── Slip4HySea00025_000.dat
	│   │       ├── Slip4HySea00026_000.dat
	│   │       ├── Slip4HySea00027_000.dat
	│   │       ├── Slip4HySea00028_000.dat
	│   │       ├── Slip4HySea00029_000.dat
	│   │       ├── Slip4HySea00030_000.dat
	│   │       ├── Slip4HySea00031_000.dat
	.............................................
	    └── 7_5000
	        └── WellsCopp1994_Normal
	            ├── output_file.log
	            ├── Slip4HySea00005_001.dat
	            ├── Slip4HySea00005_001.json
	            ├── Slip4HySea00005_002.dat
	            ├── Slip4HySea00005_002.json
	            ├── Slip4HySea00005_003.dat
	            ├── Slip4HySea00005_003.json
	            ├── Slip4HySea00005_004.dat
	            ├── Slip4HySea00005_004.json
	            ├── Slip4HySea00005_005.dat
	            ├── Slip4HySea00005_005.json
	            ├── Slip4HySea00007_001.dat
	            ├── Slip4HySea00007_001.json
	            ├── Slip4HySea00007_002.dat
	            ├── Slip4HySea00007_002.json
	            ├── Slip4HySea00007_003.dat
	            ├── Slip4HySea00007_003.json
	            ├── Slip4HySea00007_004.dat
	            ├── Slip4HySea00007_004.json
	            ├── Slip4HySea00007_005.dat
	            ├── Slip4HySea00007_005.json
	            ├── Slip4HySea00011_001.dat
	            ├── Slip4HySea00011_001.json
	            ├── Slip4HySea00011_002.dat
	            ├── Slip4HySea00011_002.json
	            ├── Slip4HySea00011_003.dat
	            ├── Slip4HySea00011_003.json	
	
For smaller magnitudes not enough cells are defined to build stochastic slip distributions and homogeneous (or modulated by rigidity) slip distributions are computed with file names `Slip4HySeaXXXXX_000.dat`. Those files are in the standard format used as input by the software [Tsunami-HySea](https://edanya.uma.es/hysea/) which is one of the most widely used tsunami simulators within the community. For larger magnitude also geojson files are computed and the second index indicated the numbering of stochastic distributions for each rupture area. Here below an example for one of the `Slip4HySeaXXXXX_00Y.dat` file

	  LON1     LAT1    DEPTH1(km)      LON2    LAT2    DEPTH2(km)      LON3    LAT3    DEPTH3(km)      RAKE    SLIP(m)
	   12.869199   38.554981    8.857142   12.812249   38.546745    8.857142   12.850728   38.512363   11.142860  -90.000000    9.104819
	   12.812249   38.546745    8.857142   12.869199   38.554981    8.857142   12.807246   38.586212    6.571429  -90.000000    8.936686
	   12.812249   38.546745    8.857142   12.794422   38.503761   11.142860   12.850728   38.512363   11.142860  -90.000000    9.739596
	   12.907272   38.519966   11.142860   12.869199   38.554981    8.857142   12.850728   38.512363   11.142860  -90.000000    7.123791
	   12.869199   38.554981    8.857142   12.863686   38.594292    6.571429   12.807246   38.586212    6.571429  -90.000000    8.475728
	   12.750983   38.577427    6.571429   12.812249   38.546745    8.857142   12.807246   38.586212    6.571429  -90.000000    9.024770
	   12.794422   38.503761   11.142860   12.812249   38.546745    8.857142   12.755468   38.537846    8.857142  -90.000000   10.431090
	   12.794422   38.503761   11.142860   12.855949   38.472931   13.428571   12.850728   38.512363   11.142860  -90.000000    8.040796
	   12.869199   38.554981    8.857142   12.907272   38.519966   11.142860   12.926314   38.562519    8.857142  -90.000000    6.726944
	   12.913015   38.480690   13.428571   12.907272   38.519966   11.142860   12.850728   38.512363   11.142860  -90.000000    4.930675
	   12.863686   38.594292    6.571429   12.869199   38.554981    8.857142   12.926314   38.562519    8.857142  -90.000000    7.812915
	   12.863686   38.594292    6.571429   12.825430   38.629196    4.285714   12.807246   38.586212    6.571429  -90.000000    8.009076
	   12.812249   38.546745    8.857142   12.750983   38.577427    6.571429   12.755468   38.537846    8.857142  -90.000000    9.790553
	   12.768585   38.620533    4.285714   12.750983   38.577427    6.571429   12.807246   38.586212    6.571429  -90.000000    9.096627
	   12.794422   38.503761   11.142860   12.755468   38.537846    8.857142   12.738218   38.494797   11.142860  -90.000000    9.190249
	   12.794422   38.503761   11.142860   12.799115   38.464211   13.428571   12.855949   38.472931   13.428571  -90.000000    7.712521
	   12.855949   38.472931   13.428571   12.913015   38.480690   13.428571   12.850728   38.512363   11.142860  -90.000000    4.808530
	   12.907272   38.519966   11.142860   12.963851   38.527435   11.142860   12.926314   38.562519    8.857142  -90.000000    5.194019
	   12.907272   38.519966   11.142860   12.913015   38.480690   13.428571   12.970135   38.488228   13.428571  -90.000000    3.787161
	   12.863686   38.594292    6.571429   12.926314   38.562519    8.857142   12.920271   38.601757    6.571429  -90.000000    6.774910
	   12.825430   38.629196    4.285714   12.863686   38.594292    6.571429   12.882524   38.636814    4.285714  -90.000000    7.181717
	   12.825430   38.629196    4.285714   12.768585   38.620533    4.285714   12.807246   38.586212    6.571429  -90.000000    8.382802
	   12.755468   38.537846    8.857142   12.750983   38.577427    6.571429   12.694775   38.568462    6.571429  -90.000000   10.362476
   .........................................................................................................................................

## 4 Post-process

 The slip distributions can be easily plotted by simple personal scripts. while the geojson files can be uploaded to Qgis (see Figure 2) or to whatever webservice using the geojson standard (e.g. [kepler.gl/](https://kepler.gl/) 

![Screenshot of a slip distributions imported to QGis](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/utils/sz_slabs/Image_qgis2.jpg)
*Figure 2 - Screenshot of a slip distributions imported to QGis* 
 
 Beyond that, in the folder [utils](https://github.com/antonioscalaunina/pyANTI-FASc/tree/main/utils) there is the script [plot_slip_distribution.py](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/utils/plot_slip_distribution.py). It might be run (still in the *antifasc* Conda environment) with the following command:
 
	python plot_slip_distribution.py

 This script will ask which class the user wants to plot as shown in the example below:
 First select among all the outputs produced so far you folder. In the shown example the ITCF00G_Hazard_slip_ITC/ folder is the 24th

 	Current folder is '../output/'

	Choose your event directory between:
	1. Irpinia_W&C1994_M69_E1533_N4085_slip_IrT/
	2. Calabrian_test_3_M70_E1800_N3760_slip_CaL/
	3. Messina_W&C1994_M70_E1569_N3817_slip_MeS/
	4. Calabrian_test_M70_E1827_N3800_slip_CaL/
	5. Samos_test_QGis_fix_mesh_M70_E2675_N3790_slip_SaM/
	6. Samos_rough_M70_E2682_N3784_slip_SaM/
	7. Irpinia_boh_3_M70_E1550_N4060_slip_IrP/
	8. Irpinia_test_QGis_M70_E1557_N4056_slip_IrP/
	9. Calabrian_test_static_M80_E1800_N3760_slip_CaA/
	10. Irpinia_refined_M70_E1550_N4060_slip_IrR/
	11. Samos_test_QGis_M70_E2675_N3790_slip_SaM/
	12. Samos_test_QGis_1_M70_E2675_N3790_slip_SaM/
	13. Irpinia_boh_M70_E1510_N4010_slip_IrP/
	14. Irpinia_refined3_M70_E1550_N4060_slip_IrR/
	15. Irpinia_test_M70_E1550_N4060_slip_IrT/
	16. Calabrian_test_7_M70_E1800_N3760_slip_CaA/
	17. Irpinia_refined2_M70_E1550_N4060_slip_IrR/
	18. Calabrian_test_8_M80_E1800_N3760_slip_CaA/
	19. Irpinia_boh_2_M70_E1480_N4080_slip_IrP/
	20. Irpinia_refined25_M70_E1550_N4060_slip_IrR/
	21. Samos_test_QGis_fix_thin_M70_E2675_N3790_slip_SaM/
	22. Calabrian_test_2_M70_E1800_N3760_slip_CaL/
	23. Sicily_rough_M70_E1311_N3850_slip_ITC/
	24. ITCF00G_Hazard_slip_ITC/
	25. Calabrian_test_9_M80_E1800_N3760_slip_CaA/
	26. Tohoku_test_json_M90_E14237_N3832_slip_KuJ/
	27. Irpinia_refined5_M70_E1550_N4060_slip_IrR/
	28. Irpinia_boh_4_M70_E1550_N4060_slip_IrP/
	29. Samos_test_QGis_large_M70_E2675_N3790_slip_SaM/
	30. Tohoku_test_M90_E14237_N3832_slip_KuJ/
	Insert a number between 1 and 30:
	
	24

Then we select either variable_mu of homogeneous_mu case

	Current folder is '../output/ITCF00G_Hazard_slip_ITC/'

	Choose your rigidity distribution directory between:
	1. variable_mu/
	2. homogeneous_mu/
	Insert a number between 1 and 2:
	
	1

Then we select one of the magnitude bin
	
	Current folder is '../output/ITCF00G_Hazard_slip_ITC/variable_mu/'

	Choose your magnitude directory between:
	1. 5_6000/
	2. 7_5000/
	3. 5_8000/
	4. 6_8000/
	5. 7_2000/
	6. 6_7000/
	7. 6_2000/
	8. 7_0000/
	9. 6_9000/
	10. 6_3000/
	11. 6_1000/
	12. 6_4000/
	13. 5_7000/
	14. 5_5000/
	15. 7_1000/
	16. 6_5000/
	17. 7_4000/
	18. 7_3000/
	19. 6_0000/
	20. 5_9000/
	21. 6_6000/
	Insert a number between 1 and 21:
	
	2

A single scaling law is defined:

	Current folder is '../output/ITCF00G_Hazard_slip_ITC/variable_mu/7_5000/'

	There is only one scaling law directory
	
	../output/ITCF00G_Hazard_slip_ITC/variable_mu/7_5000/WellsCopp1994_Normal/
	Number of files:  160
	[#####                         ]  17%



Within the selected folder, for each Slip4Hysea*.dat file, a new file will be produced that is an interactive maps in HTML format (example [here](https://antonioscalaunina.github.io/pyANTI-FASc/utils/Slip4HySea00136_005.html)). 

As already outlined in the [README](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/README.md), this example can be run through the Jupyter Notebook available [here](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/bin/antifasc_main.ipynb). For example, you can install an IDE like [Visual Studio Code](https://code.visualstudio.com/download)  to run the Jupyter Notebook on Ubuntu, macOS, or Windows. 
**IMPORTANT**: To ensure everything works correctly, it's essential to properly install all the needed components as outlined in the README.


