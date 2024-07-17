The module named *k223d* generates for each rupture area computed in the *Rupture_areas_OF* module a prescribed number of k2 stochastic slip distributions.

This module is based on the code k223d available at *https://github.com/s-murfy/k223d*. The composite source model technique there implemented is in turn based on the slipk2 code available at: *https://github.com/andherit/slipk2*. The distributed module also contains the software computing the distances over non-regular mesh surfaces through a double-lateration scheme. The kernel of this software is distributed in the github repository *https://github.com/andherit/trilateration* and described in the paper Herrero & Murphy (2018, GJI). The input configurations are read through the use of the module *forparse* available at the github repository *https://github.com/andherit/forparse*

In this repository, the software *k223d* has been modified to:

 - Compute sets of *n* stochastic slip distributions over precomptued rupture areas, and using precomputed inter-nodes distances.
 
 - Compute slip distributions accounting for a variable rigidity across the mesh

 - Write output files in the standard format of initial conditions for tsunami wave propagation simulators such as Tsunami-HySeA (see as reference *https://edanya.uma.es/hysea/index.php/models/tsunami-hysea*)

This module also makes use:

 - the module utm_geo from the repository *specfem3d* available at: *https://github.com/geodynamics/specfem3d/blob/master/src/shared/utm_geo.f90*.

 - A Fortran version of the script tiConnect_2D.m available in the github repository nodal_dg at https://github.com/tcew/nodal-dg/tree/master/Codes1.1/Codes2D/tiConnect2D.m. This script is redistributed in this distribution along with its copyright and licence notice (see in the folder *preprocess*).

BIBILIOGRAPHY

Herrero and Murphy (2018), Self-similar slip distributions on irregular shaped faults, Geophysical Journal International, DOI: 10.1093/gji/ggy104
