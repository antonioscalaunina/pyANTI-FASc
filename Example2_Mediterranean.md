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

This file [input_Sicily.json](https://github.com/antonioscalaunina/pyANTI-FASc/tree/main/config_files/Parameters/input_Sicily.json) contained in the [config_files/Parameters](https://github.com/antonioscalaunina/pyANTI-FASc/tree/main/config_files/Parameters) folder manages the main configuration parameters to run pyANTI-FASc. The file actually used for the run must be always named **input.json**. The default available file is set to run the Tohoku example, but it is sufficient to copy this `input_Sicily.json` into `input.json` to run this example

Here below the important settings to be managed by the user are shown. **Look carefully to the comments besides the parameters**. The parameters not shown in this example might be left unmodified. Their use and functionality will be fixed in next releases and better described in [Wiki Documentation](https://github.com/antonioscalaunina/pyANTI-FASc/wiki) currently under construction.

	{"zone_name": "ITCF00G",      # Name of the precomputed mesh to be used. The mesh downloaded from the ESFM20 service must be saved in geojson file as within the `utils/sz_slabs/' folder as shown [here](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/utils/sz_slabs/ITCF00G_mesh.json)
	"Merc_zone": 33,                  # Mercator zone for the selected slab. See the slab database and use the correct Mercator zone 
    "acronym": "ITC",		  # 3 digit acronym that is used for that slab. It can be arbitrarily chosen by the user (but must have 3 digits!). You might find suggestions into the slab database.
	"mesh_gen": 1,           # This option set to 1 means that a geojson file containing the mesh with name consistent with the "zone_name" (zone_name_mesh.json) is expected to be found in the utils/sz_slabs (see file [ITCF00G_mesh.json](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/utils/sz_slabs/ITCF00G_mesh.json))
    # If geojson file is not found or an incorrect name is provided, the run will be stopped and the following error message will inform the user 
	# ERROR: Mesh in GeoJSON format does not exist! Please check option in input.json and zone/file names
	"rake": -90,             # This option allows the user to set a value for the rake. If not provided a standard rake=90° is assigned to the whole mesh.
	# This field might alternatively contain the path to a csv file where the field rake is assigned, in such a case a number of rake angle consistent with the cell number of the mesh must be given. If not a standard rake=90° is assigned to the whole mesh
	
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

 The magnitude bins and the rupture geometries (defined by the selected scaling laws) are set in the input file [scaling_relationship.json](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/config_files/Parameters/scaling_relationship.json) contained in the [config_files/Parameters](https://github.com/antonioscalaunina/pyANTI-FASc/tree/main/config_files/Parameters) folder. 
The file actually used for the run must be always named **scaling_relationship.json**. The default available file is set to run this example, for which we use a selection similar to the one proposed in the framework of the project TSUMAPS-NEAM (see Basili et al. 2021) using the Strasser et al. (2010) and the Murotani et al.(2013) scaling relationships. However, it can be easily modified to run the software with different magnitude binnings and scaling laws. In the example below the structure of this file, **look carefully at the comments beside to properly set the values**:

    { 
    "Magnitude_bins": {                                     # Within this section the number of magnitude bins and the magnitude bins are defined
    "number_bins" : 32, 


