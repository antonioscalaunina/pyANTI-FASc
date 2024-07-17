# pyANTI-FASc
Python version of ANTI-FASc

ANTI-FASc, acronym for Automatic Numerical Tsunami Initial conditions: on-the-Fly rupture Areas and earthquake Scenarios, is a software enabling the fast computation of large ensembles of slip distributions on complex non-planar fault interfaces (Maesano et al. 2017, Sci. Rep.; Tonini et al. 2021, GJI) such as the subducting plates. These slip models can be promptly used as initial conditions for the computation of tsunami scenarios in the framework of both Seismic-Probabilistic Tsunami Harzard Assessment (see Scala et al. 2020 PAGEOPH - Basili et al. 2021 Frontiers) and for real-time Probabilisitic Tsunami Forecasting (see Selva et al. 2021 Nature). The newest release (version v1.1.0) is also available at the following DOI: *https://zenodo.org/doi/10.5281/zenodo.7101459* along with all the previous realeses. IMPORTANT: Please refer and cite the DOI, in your pubblications if you use the software for your research studies.

A wiki documentation (currently under construction) is available at the following link https://github.com/antonioscalaunina/ANTI-FASc/wiki

The software is composed by three modules + a simple postprocessing tool

1- Preprocess module including:
    
   - A mesh generator 
    
   - The computation of a matrix containing all the interdistances among nodes for the computed mesh
    

**More details about the preprocess modules can be found in the file *preprocess/README.md* and in the examples presented in the main folder**.

2 - Rupture areas computation:
    
   This part has two different use mode:
         
   - Hazard: it computes all the possible different rupture areas (according to the previous barycenters selection) in the prescribed magnitude bins
         
   - PTF: it computes all the scenarios “compatible” with estimation and uncertainty of magnitude and location for a given earthquake

**More details about the preprocess modules can be found in the file *bin/README.md* and in the examples presented in the main folder**

3 - k223d k-square slip distributions to 3D fault planes (Herrero & Murphy 2018, GJI):

   Computation of ensembles of stochastic k-square slip distributions for all the previously selected areas also accounting for other conditions (e.g. homogeneous or variable rigidity, surface slip amplification)
  
**More details about the k223d module and its original sources can be found in the file *src/k223d/README.md* (*https://github.com/antonioscalaunina/ANTI-FASc/blob/main/src/k223d/README.md*). The use of this module is shown in the examples in the main folder**.


Finally, the repository contains a simple postprocessing tool allowing to visualize some of the computed distributions through *.png* plots. To run it:

    cd utils
    
and:

    matlab -nodisplay -nosplash -nodesktop -r "run('slip_distribution_plot_AGI.m'); exit;"  #if you have a licensed version of MATLAB
    
or:

    ./run_slip_distribution_plot_AGI.sh /usr/local/MATLAB/MATLAB_Runtime/v99/     #if you have installe MATLAB Runtime
    
This script will ask in which folder you want to plot the computed distributions and save the plot files in the same output folder. More details in the examples.


Along with the codes it has been provided a large dataset of precomputed meshes in the folder *utils/sz_slabs*. These meshes has been computed either from the geometrical modelling provided by the project Slab 2.0 (see available geometry at the web-page https://www.sciencebase.gov/catalog/item/5aa1b00ee4b0b1c392e86467) or from the modelling proposed by the Geoscience Australia (available at: *https://github.com/GeoscienceAustralia/ptha*). The Mediterranean slab models are the same used in the framework of TSUMAPS-NEAM project (see *http://www.tsumaps-neam.eu/*, see also Maesano et al. 2017, Basili et al. 2021)

ACKNOWLEDGEMENTS

A special thanks to Stefano Lorito, Fabrizio Romano, Manuela Volpe, Hafize Basak Bayraktar, Jacopo Selva, Gaetano Festa and Antonio Giovanni Iaccarino for participating at the different phases of conceiving, revising, developing and testing of the current version of the platform.

Thanks to Roberto Basili, Francesco Emanuele Maesano, Mara Monica Tiberti and Gareth Davies for their precious contribution in providing detailed geometrical models for some of the slabs included in the database distributed along with this platform.

Thanks to Shane Murphy and Andre Herrero for their valuable contribution in developing and training to use the softwares composing the third module of the platform  



BIBLIOGRAPHY

Basili R. et al. (2021), The Making of the NEAM Tsunami Hazard Model 2018 (NEAMTHM18), Frontiers in Earth Science, DOI: 10.3389/feart.2020.616594 

Herrero and Murphy (2018), 	Self-similar slip distributions on irregular shaped faults, Geophysical Journal International, DOI: 10.1093/gji/ggy104

Maesano, F.E., Tiberti, M.M. and Basili, R. (2017), The Calabrian Arc: Three-dimensional modelling of the subduction interface, Scientific Reports DOI: 10.1038/s41598-017-09074-8

Murotani, S., Satake, K., and Fujii, Y. (2013), Scaling relations of seismic moment, rupture area, average slip, and asperity size for M~9 subduction-zone earthquakes, Geophysical Research Letters, 40, DOI: 10.1002/grl.50976.

Scala A. et al. (2020), Effect of Shallow Slip Amplification Uncertainty on Probabilistic Tsunami Hazard Analysis in Subduction Zones: Use of Long-Term Balanced Stochastic Slip Models, Pure and Applied Geophysics, DOI: 10.1007/s00024-019-02260-x

Selva, J., Lorito, S., Volpe, M. et al. (2021). Probabilistic tsunami forecasting for early warning. Nat Commun 12, 5677 (2021). DOI: 10.1038/s41467-021-25815-w

Strasser, F. O., Arango, M. C., & Bommer, J. J. (2010), Scaling of the source dimensions of interface and intraslab subduction-zone earthquakes with moment magnitude. Seismological Research Letters, DOI: 10.1785/gssrl.81.6.941.

Tonini et al. (2020), Importance of earthquake rupture geometry on tsunami modelling: The Calabrian Arc subduction interface (Italy) case study, Geophysical Journal International, DOI: 10.1093/gji/ggaa409
