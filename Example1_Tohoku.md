# EXAMPLE 1 Tohoku earthquake

In this brief guide, we a practical example for running ANTI-FASc. 
In this example we generate a set of scenarios based on location and magnitude of 2011, March the 11th, Tohoku earthquake with a magnitude Mw=9.1


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

######### UNDER CONSTRUCTION   ###############

   
