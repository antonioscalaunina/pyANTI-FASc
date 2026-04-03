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
	"name" : ["WC1994_Normal"],        # Names of scaling laws (must be consistent with the set number in the parameter "Scaling_law"
	"Area":  [43.6516,   52.7230,   63.6796,   76.9130,   92.8966,  112.2018,  135.5189,  163.6817,  
	197.6970,  238.7811,  288.4032,  348.3373,  420.7266,  508.1594,  	613.7620,  741.3102,  
	895.3648, 1081.4340, 1306.1709, 1577.6113, 1905.4607],  #Values of the area. They must be "number_bins" * "number" (of Scaling law). 
	"Length": [7.4131,    8.3176,    9.3325,   10.4713,   11.7490,   13.1826,   14.7911,   16.5959,  
	18.6209,   20.8930,   23.4423,   26.3027,   29.5121,   33.1131,   37.1535,   41.6869,   46.7735,  
	52.4807,   58.8844,   66.0693,   74.1310] #Values of the length. They must be "number_bins" * "number" (of Scaling law).
	}
	}


