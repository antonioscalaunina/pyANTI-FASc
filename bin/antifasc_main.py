#!/usr/bin/env python
# coding: utf-8

# In[1]:


#import libraries to use
import time
start_time = time.time()
import warnings
warnings.filterwarnings("ignore")
import argparse

import utils.slab as slab
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature


# ### 1. CREATE A SLAB OBJECT
# 
# A slab object is built with the attributes contained in the `input.json` file and  `scaling_relationship.json` file, and other attributes are computed internally based on the input information.


#Specify input json file path
#input_file='../config_files/Parameters/input.json'

parser = argparse.ArgumentParser(
    description="Run ANTIFASc from command line."
)

parser.add_argument(
    "--input",
    default="input.json",
    help=(
        "Input configuration file. "
        "If only a filename is provided, it is searched in config_files/Parameters. "
        "Accepted formats: .json, .yaml, .yml."
    )
)

args = parser.parse_args()

input_file = slab.prepare_input_file(args.input)

#Initialize an instance of the class Slab with the input files
Slab_obj=slab.Slab(input_file)


# In[4]:
Slab_obj.compute_matrix_distance()
Slab_obj.compute_area()
Slab_obj.Element2Element()

# In[5]:
#Select active barycenters for the whole subduction 
Slab_obj.active_barycenters()
#Selection of barycenters for case-study or hazard
Slab_obj.select_barycenter2()


# We can check the Mw and scaling relationships for which rupturing barycenters were selected
# In[6]:
print(f'Mw: {Slab_obj.get_magnitudes()}')
print(f'Scaling names: {Slab_obj.Name_scaling}')


# Compute the rupturing areas for every barycenter in each combination of magnitude and scaling relationship by applying the method `rupture_areas()`
# In[8]:
Slab_obj.rupture_areas()

# Likewise, we can check the number of rupture areas computed for each combination of magnitude and scaling relationship name
# In[9]:
Slab_obj.get_RuptAreas_number()


# Finally, write the output files
# In[11]:
#Write Output
Slab_obj.write_output_rupture_areas()
#generate folder tree for outputs
Slab_obj.generate_foldertree_slip()


# ### SLIP DISTRIBUTION COMPUTATION ###
# In[12]:
#compute slipe distributions
Slab_obj.slip_distribution()


# In[13]:
end_time = time.time()
elapsed = end_time - start_time
minutes = int(elapsed // 60)
seconds = elapsed % 60
print(f"Total running time: {minutes} min {seconds:.2f} sec")

