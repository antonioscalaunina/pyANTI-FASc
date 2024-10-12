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

Here below the important settings to be managed by the user are shown. Look carefully to the comments besides the parameters:

	{"zone_name": "kurilsjapan", 
	"Merc_zone": 54,
    	"acronym": "KuJ",
	"mesh_gen": 0,
   	"slab_file": "kur_slab2_dep_02.24.18.xyz",
    	"seismog_depth": 60,
    	"depth_interpolator": "v4",
    	"mesh_convex": 0.5,
    	"element_size": 12.5e3,
       	
	"Event": {
	"Name": "Tohoku_test",
	"Hypo_LonLat" : [142.369, 38.322],
	"Magnitude" : 9.0
	},
	"Configure": {
	"application": "PTF",
	"shape": "Rectangle",
	"numb_stoch": 5,
	"variable_mu": 1,
	"coupling_shallow_limit":2.5,
	"coupling_deep_limit":40.0,
	"mesh_sub_boundary": 0,
	"preprocess": 1,
	"file_baryc": 0,
	"file_baryc_name": "ScenarioProb_nsig2_Mw83_2015_0916_illapel_World.mat",
	"Magnitude_lb": 0.15,
	"Magnitude_ub": 0.15,
	"minimum_bnd_distance": 0.25,
	"minimum_interdistance": 0.1,
	"hypo_baryc_distance": 1.0,
	"Fact_area_scaling": 1,
	"Rigidity_file_logic": 0,
	"Rigidity_file": "Rigidity_variation.txt",
	"Stress_drop_var": 0,
	"Fact_rigidity": 0.5
	}
	}

######### UNDER CONSTRUCTION   ###############

   
