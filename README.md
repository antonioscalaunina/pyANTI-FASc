# pyANTI-FASc

## 1. General description

Python version of ANTI-FASc

pyANTI-FASc, acronym for python Automatic Numerical Tsunami Initial conditions: on-the-Fly rupture Areas and earthquake Scenarios, is a software allowing the fast computation of large ensembles of slip distributions on complex non-planar fault interfaces (Herrero & Murphy 2018 GJI, Maesano et al. 2017, Sci. Rep.; Tonini et al. 2021, GJI) such as the subducting plates. These slip models can be promptly used as initial conditions for the computation of tsunami scenarios in the framework of both Seismic-Probabilistic Tsunami Harzard Assessment (see Scala et al. 2020 PAGEOPH - Basili et al. 2021 Frontiers) and for real-time Probabilisitic Tsunami Forecasting (see Selva et al. 2021 Nature). The newest release (version v1.0.0) is also available at the following DOI: *https://zenodo.org/doi/10.5281/zenodo.13614657* along with all the previous realeses. IMPORTANT: Please refer the repository and cite the DOI as indicated at the zenodo webpage, in your pubblications, if you use the software for your research studies
Antonio Scala, Manuel Mojica, & rissclab-tester. (2024). antonioscalaunina/pyANTI-FASc: pyANTIFASc_1.0.0 (v1.0.0). Zenodo. *https://doi.org/10.5281/zenodo.13614657*

A wiki documentation (currently under construction) is available at the following link *https://github.com/antonioscalaunina/pyANTI-FASc/wiki*

The software is composed by three modules managed by a single standalone python executable that runs the module sequentially, managing the transfer of the required files among the three modules
A simple postprocessing tool is also provided to convert the output files in a standard georeference format and plot the slip distributions.

The three modules can be summarized as follows

1- **Preprocess module** 

This module includes:
    
   - A mesh generator that generates input mesh file from a nodes and faces discretization. A set of pre-computed discretization of the main subducting slabs worldwide is also provided (based on Slab 2.0 see available geometries at the web-page https://www.sciencebase.gov/catalog/item/5aa1b00ee4b0b1c392e86467) 
    
   - A preliminary computation of inter-distance between nodes to be used in the rupture areas computation
    

2 - **Rupture areas computation**:
    
This module computes a set of possible rupture geometry on the selected fault mesh. It has two different use mode:
         
   - **Hazard**: it computes a large set of possible different rupture areas in all the prescribed magnitude bins to cover in a homogeneous way the whole provided meshed zone
         
   - **PTF**: it computes all the scenarios “compatible” with estimation and uncertainty of magnitude and location for a given earthquake


3 - **k223d:** 

This module, based on the original software presented in Herrero & Murphy (2018, GJI) and available at *https://github.com/s-murfy/k223d* performs:

   - A refined computation of inter-distance between nodes to be used in k-square slip distribution computation. It is based on the lateration algorithm presented in Herrero & Murphy (2018, GJI)
   
   - The computation of ensembles of stochastic k-square slip distributions for all the previously selected areas also accounting for other conditions (e.g. homogeneous or variable rigidity, surface slip amplification, variable stress-drop etc.)
  
**IMPORTANT:** More details about the k223d module and its original sources can be found in the repository of ANTI-FASc at the following link *https://github.com/antonioscalaunina/ANTI-FASc/blob/main/src/k223d/README.md*. The use of this module is shown in the examples in the main folder.

Below are the instructions for installing the software dependencies. Please refer to the examples in the main folder and to the wiki documentation **(BOTH UNDER CONSTRUCTION)** for further details regarding the code functionality, the configuration of input files and the database of precomputed mesh discretizations

## 2 Installation

### 2.1 Linux & Windows WSL

### 2.2 Windows through conda GUI

If you want to run the code in a fully Windows environment, e.g. using a conda GUI, like Anaconda Navigator, you should follow these steps as outlined:

#### 2.2.1 Create a conda environment

1 - Download the repository at the main page *https://github.com/antonioscalaunina/pyANTI-FASc/tree/main* or with the direct link *https://github.com/antonioscalaunina/pyANTI-FASc/archive/refs/heads/main.zip* and unzip it.
2 - Open a PowerShell, within the Conda GUI you are using, and within the main folder of the repository type the following command:

    conda env create -f ANTIFASc_Win.yml

this command will create the conda environment antifasc (installing Python 3.9.16) and will install all the needed libraries and dependencies within it. Then enter into the environment typing

    conda activate antifasc

or searching for the environment antifasc through the menu *Environments* of the GUI.


## 3 ACKNOWLEDGEMENTS

A special thanks to Stefano Lorito, Fabrizio Romano, Manuela Volpe, Hafize Basak Bayraktar, Jacopo Selva, Gaetano Festa and Antonio Giovanni Iaccarino for participating at the different phases of conceiving, revising, developing and testing of the current version of the platform.

Thanks to Roberto Basili, Francesco Emanuele Maesano, Mara Monica Tiberti and Gareth Davies for their precious contribution in providing detailed geometrical models for some of the slabs included in the database distributed along with this platform.

Thanks to Shane Murphy and Andre Herrero for their valuable contribution in developing and training to use the softwares composing the third module of the platform  





## 4 BIBLIOGRAPHY

Basili R. et al. (2021), The Making of the NEAM Tsunami Hazard Model 2018 (NEAMTHM18), Frontiers in Earth Science, DOI: 10.3389/feart.2020.616594 

Herrero and Murphy (2018), 	Self-similar slip distributions on irregular shaped faults, Geophysical Journal International, DOI: 10.1093/gji/ggy104

Maesano, F.E., Tiberti, M.M. and Basili, R. (2017), The Calabrian Arc: Three-dimensional modelling of the subduction interface, Scientific Reports DOI: 10.1038/s41598-017-09074-8

Murotani, S., Satake, K., and Fujii, Y. (2013), Scaling relations of seismic moment, rupture area, average slip, and asperity size for M~9 subduction-zone earthquakes, Geophysical Research Letters, 40, DOI: 10.1002/grl.50976.

Scala A. et al. (2020), Effect of Shallow Slip Amplification Uncertainty on Probabilistic Tsunami Hazard Analysis in Subduction Zones: Use of Long-Term Balanced Stochastic Slip Models, Pure and Applied Geophysics, DOI: 10.1007/s00024-019-02260-x

Selva, J., Lorito, S., Volpe, M. et al. (2021). Probabilistic tsunami forecasting for early warning. Nat Commun 12, 5677 (2021). DOI: 10.1038/s41467-021-25815-w

Strasser, F. O., Arango, M. C., & Bommer, J. J. (2010), Scaling of the source dimensions of interface and intraslab subduction-zone earthquakes with moment magnitude. Seismological Research Letters, DOI: 10.1785/gssrl.81.6.941.

Tonini et al. (2020), Importance of earthquake rupture geometry on tsunami modelling: The Calabrian Arc subduction interface (Italy) case study, Geophysical Journal International, DOI: 10.1093/gji/ggaa409
