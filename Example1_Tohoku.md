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


######### UNDER CONSTRUCTION   ###############

   
