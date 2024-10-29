"""
June 2024

"""

###############################################################
#Import libraries to be used
import os
import sys
import glob
import json
import shutil
import subprocess
import platform
import numpy as np
from tqdm import tqdm
import scipy
from scipy.linalg import lstsq
from scipy.spatial import Delaunay
from scipy.io import loadmat
import cartopy.crs as ccrs

import utm
import alphashape
from shapely.geometry import Point, Polygon
from shapely.vectorized import contains
from geopy.distance import geodesic
import geopy.distance as distance

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.tri as tri


def find_main_dir():
    current_dir = os.getcwd()
    while True:
        # Check if the current directory contains the 'bin' folder
        if os.path.exists(os.path.join(current_dir, 'bin', 'utils')):
            return current_dir
        # Move one level up in the directory structure
        new_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
        if new_dir == current_dir:
            # We have reached the root directory and did not find 'main'
            raise RuntimeError("Could not find 'main' directory, please run your script inside ANTIFASc")
        current_dir = new_dir

main_dir = find_main_dir()
utils_dir= os.path.join(main_dir, 'bin', 'utils')
sys.path.append(utils_dir)

import mesh as mesh
import Rupture_areas_utils as rupture
import plot_utils as plotutils



###############################################################



class Slab:
    """
    
    """


    def __init__(self, config_file, scaling_file):
        """
        Initialize the Slab object by reading configuration and scaling files, and the mesh.
        """
        
        
        print('reading '+ config_file.split('/')[-1] + ' file')
        read_inputfile(self, config_file)
        print('reading '+ scaling_file.split('/')[-1] + ' file')
        read_scalingrel_file(self, scaling_file)

        if self.Preprocess_logic:
            mesh.faces_nodes_2_mesh_file(config_file)

        print('reading mesh')
        read_mesh(self)




    def _check_attribute(self, attr_name):
        """
        Check if the attribute attr_name exists and raise an error if it doesn't.

        Parameters:
        attr_name (str): Name of the attribute to check.
        """
        if not hasattr(self, attr_name):
            if attr_name == 'barycenter':
                raise AttributeError(f"The attribute '{attr_name}' has not been computed yet. Please run the method select_barycenter2() first")
            elif attr_name=='Rupturing_areas':
                raise AttributeError(f"The attribute '{attr_name}' has not been computed yet. Please run the method rupture_areas() first")
            else:
                raise AttributeError(f"The attribute '{attr_name}' has not been computed yet. Please ensure it is computed before accessing it.")
        



    def active_barycenters(self):
        """
        Get indices of barycenters located at greater distances than min_bnd_dist*WidthSL from the slab boundary
        
        index_active (): Indices of active barycenters for the whole subduction
        """
        
        #Compute distances between barycenters and the slab boundary
        p0, distanceJB=rupture.compute_distance2fault2D(self.bnd_mesh.copy(),self.barycenters_all.copy(),self.Sub_boundary_logic,self.Merc_zone,self.hemisphere)
        # Initialize index_active as a list of lists
        index_active = [[None for _ in range(self.N_scaling)] for _ in range(len(self.Magnitude))]
        # Fill index_active with indices of barycenters separated enough from the boundary
        for i in range(len(self.Magnitude)):
            for j in range(self.N_scaling):  # 1 Murotani - 2 Strasser
                index_active[i][j] = np.where(distanceJB > self.min_bnd_dist * self.WidthSL[i, j])[0]
        self.index_active = index_active


    
    def select_barycenter2(self):
        """
        selection of barycenters for case-study or hazard
        """
        self._check_attribute('index_active')
        
        #assuming there is not preprocess part anymore
        index_active = self.index_active
        barycenter = [[[] for _ in range(self.N_scaling)] for _ in range(len(self.index_magnitude))]
        
        print('Barycenter selection')
        
        for l, i in enumerate(self.index_magnitude):

            #if self.application == 'PTF':
            print(f'Magnitude bin # {i} - Mw={self.Magnitude[i]:.4f}')
        
            for j in range(self.N_scaling):
                 if len(index_active[i][j])==0 or self.int_dist * self.LengthSL[i, j] < 0.5 * self.elem_size / 1e3:
                    barycenter[l][j] = index_active[i][j] 
                 else:
                    barycenter[l][j] = []
                    selected_points = []
        
                    for k in index_active[i][j]:
                        point_LL = self.barycenters_all[k]
                        point=np.zeros([3])
                        point[0],point[1]=mesh.from_latlon(point_LL[1], point_LL[0], zone=self.Merc_zone)
                        #point[0],point[1]=mesh.ll2xy(point_LL[1], point_LL[0], zone=self.Merc_zone,hemisphere=self.hemisphere)
                        point[2]=point_LL[2]
        
                        if not selected_points:
                            selected_points.append(point)
                            barycenter[l][j].append(k)
                        else:
                            distances = np.sqrt(np.sum((np.array(selected_points) - np.array(point))**2, axis=1))
                            if np.min(distances) >= self.int_dist * self.LengthSL[i, j] * 1e3:
                                selected_points.append(point)
                                barycenter[l][j].append(k)
        
        
        ind_aux=barycenter.copy()
        barycenter = [[[] for _ in range(self.N_scaling)] for _ in range(len(self.index_magnitude))]
        if self.application == 'PTF':  # Discrimination is the distance from the hypocenter
            hypo_GEO = self.hypo_GEO
            print('Barycenter selection (PTF)')
            for l, i in enumerate(self.index_magnitude):
                print(f'Magnitude bin # {i} - Mw={self.Magnitude[i]:.4f}')
                #l += 1
                for j in range(self.N_scaling):
                    #kk = 0
                    for k in ind_aux[l][j]:
                        # Calculate the distance between the barycenter and the hypocenter
                        distance = geodesic((self.barycenters_all[k, 1], self.barycenters_all[k, 0]), (self.hypo_GEO[1], self.hypo_GEO[0])).meters
                        if distance < self.hypo_baryc_dist * self.LengthSL[i, j] * 1e3:
                        
                            barycenter[l][j].append(k)
                            #kk += 1
                    #if kk == 1:
                        #barycenter[l-1][j] = []

        elif self.application == 'Hazard':  # Discrimination is the distance from the active barycenter
            barycenter = ind_aux
        
        self.barycenter = barycenter
        self.ind_aux=ind_aux

   
    

    def Element2Element(self):
        """
        Compute or load element-to-element connectivity.
        """
        """if self.Preprocess_logic:
            file_EToE='../config_files/Connection_cell/EToE_'+self.zone_code+'.mat'
            data_EToE = loadmat(file_EToE)
            self.EToE=data_EToE['EToE']
        else:"""
        EToE,EToF=mesh.compute_connectivity(self.cells)
        self.EToE=EToE
        self.EToF=EToF



    def compute_matrix_distance(self):

        #if self.Preprocess_logic:
            #Matrix_distance_file='../config_files/Matrix_distances/'+self.zone_code+'_matrix_distance.bin'
            #Matrix_distance = np.fromfile(Matrix_distance_file, dtype=np.float32)
            #nodes_length = len(self.nodes)
            #Matrix_distance_shape = (nodes_length, nodes_length)
            #Matrix_distance = Matrix_distance.reshape(Matrix_distance_shape)
            #self.Matrix_distance=Matrix_distance
        #else:
        self.Matrix_distance=mesh.matrix_distance_nolat2(self.nodes,self.zone_code)
    


    def compute_area(self):
        """
        Compute the total and individual cell areas.
        """
        self._check_attribute('Matrix_distance')
        Area_tot, Area_cells=mesh.compute_area(self.cells,self.Matrix_distance)
        self.Area_cells=Area_cells*1e-6



    def rupture_areas(self):
        
        self._check_attribute('barycenter')
        self._check_attribute('Matrix_distance')
        self._check_attribute('EToE')

        print('Rupturing area computation')
        if self.application == 'Hazard':
            n=len(self.index_magnitude)
            progress_bar = tqdm(total=n, desc="Computing Rupturing areas")

        Rupturing_areas = [[[] for _ in range(self.N_scaling)] for _ in range(len(self.index_magnitude))]
        
        for l, i in enumerate(self.index_magnitude):

            if self.application == 'PTF':
                print(f'Magnitude bin # {i} - Mw={self.Magnitude[i]:.4f}')
            elif self.application=='Hazard':
                progress_bar.update(1)

        
            for j in range(self.N_scaling):

                if len(self.barycenter[l][j])!=0:
                    Rupturing_areas[l][j]=rupture.Rupture_area_cells(self,l,j)
                
        Slab.Rupturing_areas=Rupturing_areas

        print('Rupturing areas computed!')

    def write_output_rupture_areas(self):
        
        self._check_attribute('Rupturing_areas')
        self._check_attribute('Area_cells')
        
        #create output folders 
        os.makedirs(os.path.join(main_dir,self.namefolder))
        os.makedirs(os.path.join(main_dir,self.namefolder_slip))
        write_output(self.Rupturing_areas, self.Magnitude, self.index_magnitude, self.SPDF_all, self.Name_scaling, self.Area_cells, self.namefolder)
   


    def generate_foldertree_slip(self):
        generate_foldertree_slip(self.namefolder_slip, self.index_magnitude, self.Magnitude, self.Name_scaling)


    def get_nodes_coords(self,barycenter_indices=[]):
        """
        Get nodes coordinates of the mesh
        """
        if len(barycenter_indices)!=0:
            lon= self.nodes[barycenter_indices,0].copy()
            lat=self.nodes[barycenter_indices,1].copy()
            depth=self.nodes[barycenter_indices,2].copy()
        else:
            lon=self.nodes[:,0].copy()
            lat=self.nodes[:,1].copy()
            depth=self.nodes[:,2].copy()
    
        return lon, lat, depth

    def get_elements(self, barycenter_indices=[]):
        """
        Get indices of the vertices of the mesh's elements
        """
        
        if len(barycenter_indices)!=0:
            cells=self.cells[barycenter_indices,:]
        else:
            cells=self.cells
        
        return cells-1 #to set indices in Pyhton

    def get_barycenter_coordinates(self,barycenter_indices=[]):
        
        if len(barycenter_indices)!=0:
            lon=self.barycenters_all[barycenter_indices,0]
            lat=self.barycenters_all[barycenter_indices,1]
            depth=self.barycenters_all[barycenter_indices,2]
        else:
            lon=self.barycenters_all[:,0]
            lat=self.barycenters_all[:,1]
            depth=self.barycenters_all[:,2]
    
        return lon, lat, depth




    def get_magnitudes(self):
        """
        Return magnitudes that are actually used
        """
        return np.array(self.Magnitude)[self.index_magnitude]




    def get_magnitude_index(self, magnitude):
        if magnitude is None:
            return None
        try:
            return np.where(np.array(self.get_magnitudes()) == magnitude)[0][0]
        except IndexError:
            print(f"Magnitude {magnitude} not found.")
            print(f'Please choose one of the following magnitudes: {self.get_magnitudes()}')
            return None
    



    def get_scaling_index(self, Name_scaling):
        if Name_scaling is None:
            return None
        try:
            return np.where(np.array(self.Name_scaling) == Name_scaling)[0][0]
        except IndexError:
            print(f"Scaling Name {Name_scaling} not found.")
            print(f'Please choose one of the following names: {self.Name_scaling}')  
            return None



    
    def get_indices_scalingrel(self, magnitude=None, Name_scaling=None):
        if magnitude is None and Name_scaling is None:
            print('Please provide magnitude and/or Name_scaling')
            
    
        mag_idx = self.get_magnitude_index(magnitude)
        scaling_idx = self.get_scaling_index(Name_scaling)
    
        return mag_idx, scaling_idx
    
    


    def get_RuptAreas_number(self, magnitude=None, Name_scaling=None):
        """
        Get the number of rupturing areas for the specified magnitude and scaling relationship.
        
        Parameters:
        magnitude (float, optional): Magnitude value.
        Name_scaling (str, optional): Name of the scaling relationship.
        
        Returns:
        text: Number of rupturing areas for each required configuration of magnitude and scaling rel.
        """
        self._check_attribute('Rupturing_areas')
        
        if magnitude is None and Name_scaling is None:
            # Print all configurations
            for i, mag in enumerate(self.get_magnitudes()):
                for j, scaling in enumerate(self.Name_scaling):

                    event=self.Rupturing_areas[i][j]
                    N = sum([e.get('true', None) for e in event])
                    N_all = len(self.barycenter[i][j])
                    print(f'Mw={mag}, Name scaling: {scaling}, N={N}, N_all={N_all}')
            return
    
        mag_idx, scaling_idx = self.get_indices_scalingrel(magnitude, Name_scaling)
        if Name_scaling is None:
            for j, scaling in enumerate(self.Name_scaling):

                event=self.Rupturing_areas[mag_idx][j]
                N = sum([e.get('true', None) for e in event])
                #N = len(self.barycenter[mag_idx][j])
                print(f'Name scaling: {scaling}, N={N}')
        elif magnitude is None:
            for i, mag in enumerate(self.get_magnitudes()):
                event=self.Rupturing_areas[i][scaling_idx]
                N = sum([e.get('true', None) for e in event])
                #N = len(self.barycenter[i][scaling_idx])
                print(f'Mw={mag}, N={N}')
        else:
            event=self.Rupturing_areas[mag_idx][scaling_idx]
            N = sum([e.get('true', None) for e in event])
            #N = len(self.barycenter[mag_idx][scaling_idx])
            print(f'N={N}')

    

    def slip_distribution(self):
        """
        Run the scripts to compute the slip distributions. 
        
        """
        curr_path=os.getcwd()
        normalized_path = os.path.normpath(curr_path)
        curr_dir=os.path.basename(normalized_path)

        #if curr_dir == 'ANTI-FASc':

        if self.variable_mu==0:
            print('Computing slip distributions for the homogeneous case')
            run_homo(self)
        elif self.variable_mu==1:
            print('Computing slip distributions for the homogeneous and variable rigidity cases ')
            run_homo(self)
            run_var(self)

        #else:
           #print('Please go to the main directory to run this method')


    def plot_basemap(self, fig, subplots):
        """
        Create GeoAxes and plot basemap considering extention of the slab to plot
        """
        lon, lat, depth = self.get_nodes_coords()
        lon_0=(max(lon)+min(lon))/2 
       
        #Create GeoAxes
        ax=fig.add_subplot(subplots, projection=ccrs.PlateCarree(central_longitude=lon_0))
        #Plot basemap
        fig, ax=plotutils.set_dfbasemap(fig,ax)
        #set plot extent
        ax=plotutils.set_extent(lon,lat,ax)


        return fig, ax


    def plot_slab(self,ax,fig,colorbar=False,save=False):
        """
        Plot the slab mesh on the given axis.
        
        Parameters:
        ax (GeoAxes): Matplotlib axis to plot on.
        colorbar (bool, optional): Whether to include a colorbar. Defaults to False.
        """
        #lon,lat,depth=self.get_barycenter_coordinates()
        lon, lat, depth = self.get_nodes_coords()
        triang=self.get_elements()

        ax=plotutils.set_extent(lon,lat,ax)
        #s=plotutils.get_marker_size(self,ax,fig)
        #get central longitude and adjust longitudes
        lon_0=ax.projection.proj4_params['lon_0']    
        lon-=lon_0
    
        #plot Slab mesh
        #sc = ax.scatter(lon, lat, c=-depth/1000, cmap='viridis_r', s=s,lw=0,transform=ccrs.PlateCarree())
        tpc = ax.tripcolor(lon,lat,triang,-depth/1000,shading='gouraud',cmap='viridis_r')
    
        if colorbar:
            # Add a colorbar
            cbar=plotutils.set_cbar(tpc,ax,fig,label='Depth [km]')

        if save:
            fig.savefig(f'slab_{self.zone_code}.png')

        return ax, tpc, cbar
        
   

    def plot_SPDF(self,ax,fig,colorbar=False, save=False):
        """
        Plot the SPDF values on the given axis.
        
        Parameters:
        ax (GeoAxes): Matplotlib axis to plot on.
        colorbar (bool, optional): Whether to include a colorbar. Defaults to False.
        """
        #lon,lat,depth=self.get_barycenter_coordinates()
        lon, lat, depth = self.get_nodes_coords()
        triang=self.get_elements()
        #get central longitude and adjust longitudes
        lon_0=ax.projection.proj4_params['lon_0']
        lon-=lon_0
        #set plot extent
        #ax=plotutils.set_extent(lon,lat,ax)
        #plot SPDF
        SPDF=self.SPDF_all

        #sc = ax.scatter(lon, lat, c=SPDF, cmap='jet', s=s,lw=0, transform=ccrs.PlateCarree())
        tpc = ax.tripcolor(lon,lat,triang, SPDF,shading='flat',cmap='jet')
    
        if colorbar:
            # Add a colorbar
            cbar=plotutils.set_cbar(tpc,ax,fig,label='SPDF',scale=4)
        
        if save:
            fig.savefig(f'SPDF_{self.zone_code}.png')

        return ax, tpc, cbar
    



    def plot_barycenters_mag(self,magnitude,Name_scaling,ax,fig):
        """
        Plot the barycenters for the specified magnitude and scaling relationship on the given axis.
        
        Parameters:
        magnitude (float): Magnitude value.
        Name_scaling (str): Name of the scaling relationship.
        ax (GeoAxes): Matplotlib axis to plot on.
        """
        self._check_attribute('barycenter')
        #check if the given magnitude and name_scaling are indeed in the parameters used
        
        mag_indx,scaling_indx=self.get_indices_scalingrel(magnitude=magnitude, Name_scaling=Name_scaling)
        
        if mag_indx is not None and scaling_indx is not None:
            #find barycenter indices and coordinates for the given magnitude and name_scaling
            selected_bary_ind=self.barycenter[mag_indx][scaling_indx]
            lon, lat, depth = self.get_barycenter_coordinates(selected_bary_ind)

            lon_0=ax.projection.proj4_params['lon_0']
            lon-=lon_0
            
            #plot barycenters
            _,s=plotutils.get_marker_size(self,ax,fig)
            #ax.plot(lon,lat,'ko',markersize=s,label='Barycenters')
            ax.scatter(lon,lat,marker='o', color='black', edgecolor=None, s=10*s)
        
        if self.application=='PTF':
            sh=s*300
            self.plot_hypo(ax,sh)
    
    
    def plot_rupture_area(self, magnitude, Name_scaling, N_area, ax,fig, colorbar=False ):

        """
        Plot the specified rupture area on the given axis.
        
        Parameters:
        magnitude (float): Magnitude value.
        Name_scaling (str): Name of the scaling relationship.
        N_area (int): Index of the rupture area to plot.
        ax (GeoAxes): Matplotlib axis to plot on.
        colorbar (bool, optional): Whether to include a colorbar. Defaults to False.
        """
        self._check_attribute('Rupturing_areas')
        
        #check if the given magnitude and name_scaling are indeed in the parameters used
        mag_indx,scaling_indx=self.get_indices_scalingrel(magnitude=magnitude, Name_scaling=Name_scaling)
    
        if mag_indx is not None and scaling_indx is not None:
            
            #lon, lat, depth = self.get_barycenter_coordinates()
            #sc_sl = ax.scatter(lon, lat, c='white',s=1)
            lon, lat, depth = self.get_nodes_coords()
            triang=self.get_elements()
            lon_0=ax.projection.proj4_params['lon_0']
            lon-=lon_0

            #ax=plotutils.set_extent(lon,lat,ax)
    
            #plot Slab mesh
            tpc = ax.tripcolor(lon,lat,triang,-depth/1000,shading='gouraud',alpha=0.1,cmap='viridis_r')
            
            #Get rupturing areas 
            event=self.Rupturing_areas[mag_indx][scaling_indx]
            #remove false-flag areas
            #event = [e for e, m in zip(event, [e.get('true', None) for e in event]) if m]

            
            if N_area<len(event):
                
                bary_area=event[N_area]['cell']
                #lon_ra, lat_ra, depth_ra = self.get_barycenter_coordinates(bary_area)
        
                triang_ra=self.get_elements(bary_area)

                tpc = ax.tripcolor(lon,lat,triang_ra,-depth/1000,shading='gouraud',cmap='viridis_r')
                #sc_ra = ax.scatter(lon_ra, lat_ra, c=-depth_ra/1000, cmap='viridis_r', s=1, transform=ccrs.PlateCarree())
    
                if colorbar:
                    # Add a colorbar
                    cbar=plotutils.set_cbar(tpc,ax,fig,label='Depth [km]')
                    
            else:
                print(f'There are only {len(event)} rupturing areas for this Mw={magnitude} and Name_scaling:{Name_scaling} configuration, please use an index between 0 and {len(event)-1}')
    


    def plot_slip_dist(self, magnitude, Name_scaling, N_area, ax,fig,var=False, colorbar=False ):
        """
        Plot the specified slip distribution(s) on the given axis.
         
        Parameters:
        magnitude (float): Magnitude value.
        Name_scaling (str): Name of the scaling relationship.
        N_area (int): Index of the rupture area to plot.
        ax (GeoAxes): Matplotlib axis to plot on.
        var (bool, optional): whether to plot slip dists of variable mu
        colorbar (bool, optional): Whether to include a colorbar. Defaults to False.
        """
         
        event_out=self.namefolder_slip
        if var:
            rigidity='variable_mu'
        else:
            rigidity='homogeneous_mu'
        mw_string = f"{magnitude:6.4f}".replace('.', '_', 1)
        #get path of slip distributions
        folder_in = os.getcwd()
        folder_out = f"{main_dir}/output/{event_out}/{rigidity}/{mw_string}/{Name_scaling}/"

        mag_indx,scaling_indx=self.get_indices_scalingrel(magnitude=magnitude, Name_scaling=Name_scaling)
        
        if mag_indx is not None and scaling_indx is not None:
            
            #lon, lat, depth = self.get_barycenter_coordinates()
            #sc_sl = ax.scatter(lon, lat, c='white',s=1)
            lon, lat, depth = self.get_nodes_coords()
            triang=self.get_elements()
            lon_0=ax.projection.proj4_params['lon_0']
            lon-=lon_0

            #ax=plotutils.set_extent(lon,lat,ax)
    
            #plot Slab mesh
            tpc = ax.tripcolor(lon,lat,triang,-depth/1000,shading='gouraud',alpha=0.1,cmap='viridis_r')
    
            #Get rupturing areas 
            event=self.Rupturing_areas[mag_indx][scaling_indx]
            #remove false-flag areas
            event = [e for e, m in zip(event, [e.get('true', None) for e in event]) if m]

            
            if N_area<len(event):
                
                bary_area=event[N_area]['cell']
                triang_ra=self.get_elements(bary_area)

                ra_id=event[N_area]['cell'][0]
                slip_filename = f'{folder_out}Slip4HySea{ra_id+1:05d}_001.dat'
                slip=plotutils.get_slip(slip_filename)

                tpc = ax.tripcolor(lon,lat,triang_ra,slip,shading='flat',cmap='jet_r')
                #sc_ra = ax.scatter(lon_ra, lat_ra, c=-depth_ra/1000, cmap='viridis_r', s=1, transform=ccrs.PlateCarree())
    
                if colorbar:
                    # Add a colorbar
                    cbar=plotutils.set_cbar(tpc,ax,fig,label='Slip [m]')
                    
            else:
                print(f'There are only {len(event)} rupturing areas for this Mw={magnitude} and Name_scaling:{Name_scaling} configuration, please use an index between 0 and {len(event)-1}')
    


    def plot_hypo(self,ax,s):

        """
        Plot the hypocenter on the given axis if the application is 'PTF'.
        
        Parameters:
        ax (GeoAxes): Matplotlib axis to plot on.
        """

        if self.application=='PTF':
            lon_0=ax.projection.proj4_params['lon_0']
            ax.scatter(self.hypo_GEO[0]-lon_0,self.hypo_GEO[1], marker='*', color='yellow', edgecolor='black', s=2*s)
        else:
            print('The application is Hazard, no hypocenter is specified')
        
    def plot_boundary(self,ax):

        """
        Plot slab/fault (sub)boundary

        Parameters:
        ax (GeoAxes): Matplotlib axis to plot on.

        """
        
        lon=self.bnd_mesh[:,0].copy()
        lat=self.bnd_mesh[:,1].copy()
        lon_0=ax.projection.proj4_params['lon_0']
        lon-=lon_0
        ax.plot(lon,lat,'k-')


###----- Auxiliary Functions to build the class Slab --------###



def read_inputfile(Slab, config_file):
    # Read input from JSON file
    with open(config_file) as fid:
        Param = json.load(fid)
    
    zone_code = Param['acronym'] #Slab
    Merc_zone = Param['Merc_zone'] #Mercator projection zone
    application = Param['Configure']['application'] #Hazard: All magnitude bins - all barycenters, PTF: Magnitude and location around estimated ones
    shape = Param['Configure']['shape'] #Rectangle: rectangular shape, Circle: circular shape
    elem_size = Param['element_size'] 
    Fact_rigidity = Param['Configure']['Fact_rigidity']

    Slab.zone_code = zone_code
    Slab.Merc_zone = Merc_zone
    Slab.application = application
    Slab.shape = shape
    Slab.elem_size = elem_size
    Slab.Fact_rigidity = Fact_rigidity

    Slab.int_dist=Param['Configure']['minimum_interdistance']
    Slab.hypo_baryc_dist=Param['Configure']['hypo_baryc_distance']
    Slab.min_bnd_dist=Param['Configure']['minimum_bnd_distance']
    Slab.Fact_Area=Param['Configure']['Fact_area_scaling']
    Slab.K_SL=Param['Configure']['coupling_shallow_limit']
    Slab.K_DL=Param['Configure']['coupling_deep_limit']
    Slab.Preprocess_logic = False
    Slab.Sub_boundary_logic = False
    Slab.Stress_drop_logic = False
    Slab.Rigidity_file_logic = False
    Slab.Rigidity_file = None
    Slab.Mw = None
    Slab.hypo_GEO = None
    Slab.Baryc_file_logic=False
    Slab.Baryc_file=None
    Slab.fact_mu_z=None
    Slab.hemisphere='N'
    
    # Create Slab instance
    #Slab = Slab(zone_code, Merc_zone, application, shape, elem_size, Fact_rigidity)
    
    # Write to param_zone.dat file
    with open(os.path.join(main_dir,'param_zone.dat'), 'w') as fid:
        fid.write(f'geo zone={zone_code}\n')
        fid.write(f'mercator={Merc_zone}\n')
    
    #logical: if file with a prescribed depth variation of rigidity is given
    if Param['Configure']['Rigidity_file_logic'] == 1:
        Slab.Rigidity_file_logic = True
        Slab.Rigidity_file = Param['Configure']['Rigidity_file'] #Name of the file containing the rigidity variation with depth
    
    Slab.Preprocess_logic = Param['Configure']['preprocess'] == 1 #true: the active barycenters are already in preprocessing file
    Slab.Sub_boundary_logic = Param['Configure']['mesh_sub_boundary'] == 1 #logical: if the rupture can extend only on a part of the mesh 
    Slab.Stress_drop_logic = Param['Configure']['Stress_drop_var'] == 1 #logical: if the stress drop varies with depth
    
    # Logic for setting other attributes based on application type
    if application == 'PTF':
        Zone = Param['Event']['Name']
        Mw = Param['Event']['Magnitude']
        hypo_GEO = Param['Event']['Hypo_LonLat']
        Mw_string = f'{Mw:.1f}'.replace('.', '')
        hypo_GEO_string = [f'{abs(hypo_GEO[0]):.2f}'.replace('.', '').replace(' ', '0'),
                           f'{abs(hypo_GEO[1]):.2f}'.replace('.', '').replace(' ', '0')]
        hypo_GEO_string[0] = f"E{hypo_GEO_string[0]}" if hypo_GEO[0] >= 0 else f"W{hypo_GEO_string[0]}"
        hypo_GEO_string[1] = f"N{hypo_GEO_string[1]}" if hypo_GEO[1] >= 0 else f"S{hypo_GEO_string[1]}"
        
        Slab.Mw = Mw
        Slab.hypo_GEO = hypo_GEO
        Slab.lb_Mw = Param['Configure']['Magnitude_lb']
        Slab.ub_Mw = Param['Configure']['Magnitude_ub']
    
        namefolder = f'{Zone}_M{Mw_string}_{hypo_GEO_string[0]}_{hypo_GEO_string[1]}'
        namefolder_slip = f'{namefolder}_slip_{zone_code}'

        if Param['Configure']['file_baryc'] == 1:
            Slab.Baryc_file_logic=True
            Slab.Baryc_file=Param['Configure']['file_baryc_name']
            
    elif application == 'Hazard':
        Zone = Param['zone_name']
        namefolder = f'{Zone}_Hazard'
        namefolder_slip = f'{namefolder}_slip_{zone_code}'
    
    # Write to name_folders_file.dat file
    with open(os.path.join(main_dir,'name_folders_file.dat'), 'w') as fid:
        fid.write(f'{namefolder}\n{namefolder_slip}\n{zone_code}\n')
        fid.write(f"{Param['Configure']['numb_stoch']}\n")
        fid.write(f"{Param['Configure']['variable_mu']}\n")
    
    Slab.namefolder=namefolder
    Slab.namefolder_slip=namefolder_slip
    Slab.numb_stoch=Param['Configure']['numb_stoch']
    Slab.variable_mu=Param['Configure']['variable_mu']




def read_scalingrel_file(Slab, scaling_file):
    ## Magnitude bins and scaling laws (Strasser and Murotani size in km - km^2)

    # Load data from scaling_relationship.json
    with open(scaling_file) as fid:
        Scaling = json.load(fid)
    
    
    Magnitude = Scaling['Magnitude_bins']['Magnitude']
    N_scaling = Scaling['Scaling_law']['number']
    Name_scaling = Scaling['Scaling_law']['name']
    Area_scaling = Scaling['Scaling_law']['Area']
    Length_scaling = Scaling['Scaling_law']['Length']
    N_Magn_bins = Scaling['Magnitude_bins']['number_bins']
    
    Slab.Magnitude=Magnitude
    Slab.N_scaling=N_scaling
    Slab.Name_scaling=Name_scaling
    
    # Write to classes_scaling.dat file
    with open(os.path.join(main_dir,'config_files','Parameters','classes_scaling.dat'), 'w') as fid:
        for name in Name_scaling:
            fid.write(f'{name}\n')
    
    if len(Magnitude) != N_Magn_bins:
        print('ERROR: Number of Magnitude bin not correct')
    if len(Area_scaling) != N_Magn_bins * N_scaling:
        print('ERROR: Number of Area data not correct')
    if len(Length_scaling) != N_Magn_bins * N_scaling:
        print('ERROR: Number of Length data not correct')
    if len(Name_scaling) != N_scaling:
        print('ERROR: Number of Scaling relationship names not correct')
    
    Area_aux = np.reshape(Area_scaling, (N_scaling,N_Magn_bins))
    Length_aux = np.reshape(Length_scaling, (N_scaling,N_Magn_bins))
    Width_aux = Area_aux / Length_aux
    
    Slab.AreaSL=Area_aux.T
    Slab.WidthSL=Width_aux.T
    Slab.LengthSL=Length_aux.T 
    
    ##########################################################################
    if Slab.Stress_drop_logic:
        gamma1 = np.reshape(Scaling['Scaling_law']['gamma1'], (1, N_scaling))
        gamma2 = np.reshape(Scaling['Scaling_law']['gamma2'], (1, N_scaling))
        Slab.gamma1 = gamma1
        Slab.gamma2 = gamma2
    
    ##########################################################################
    
    # Magnitude bins 
    index_magnitude = None
    if Slab.application == 'PTF':
        
        if Slab.Baryc_file_logic:
            
            # Load from file
            """This part here has to be checked"""
            with open(os.path.join(main_dir,'config_files','PTF_selection',Slab.Baryc_file)) as fid:#which kind of file is this?
                ScenarioProb = fid.read()
            Mag_ParPS = [x[1] for x in ScenarioProb['ParScenPS']]
            index_PS = [i for i, x in enumerate(Mag_ParPS)]
            index_magnitude = [Magnitude.index(x) for x in Mag_ParPS]
        else:
            #Selection on the magnitude boundaries HERE
            lb_Mw = Slab.Mw - Slab.lb_Mw
            ub_Mw = Slab.Mw + Slab.ub_Mw
            index_magnitude = [i for i, x in enumerate(Magnitude) if lb_Mw <= x <= ub_Mw]
            
    elif Slab.application == 'Hazard':
        index_magnitude = list(range(len(Magnitude)))
    
    # Convert index_magnitude to a numpy array if needed
    if index_magnitude is not None:
        index_magnitude = np.array(index_magnitude)
    Slab.index_magnitude=index_magnitude 




### read mesh, cell barycenters and boundary of seismogenic zone and rigidity yes/no ###
def read_mesh(Slab):
    
    name_filemesh = os.path.join(main_dir,'config_files','Mesh',f'{Slab.zone_code}_mesh_15km.inp')
    with open(name_filemesh) as fid: #mesh file to be checked
        nodes, cells, _,_,_ = mesh.read_mesh_file(fid)
    
    barycenters_all = mesh.find_barycenters(nodes, cells) #compute barycenter coordinates of each cell 
    hemisphere = mesh.get_hemisphere(nodes[:,1]) #get hemisphere for further computations

    if Slab.Sub_boundary_logic:
        name_bnd = os.path.join(main_dir,'config_files','Mesh',f"{Slab.zone_code}_boundary.txt")#f"../config_files/Mesh/{Slab.zone_code}_boundary.txt"
        bnd_mesh = np.genfromtxt(name_bnd, skip_header=1)
        bnd_mesh[bnd_mesh[:, 0] < 0, 0] += 360
    else:
        nodes_plus = nodes.copy()
        nodes_plus[nodes[:, 0] < 0, 0] += 360
        hull = alphashape.alphashape(nodes[:,[0,1]],alpha=0.3)
        boundary_points = np.array(hull.exterior.coords)  # Exclude the repeated last point
        #boundary_set = set(map(tuple, boundary_points))
        #bnd = [i for i, point in enumerate(nodes[:,[0,1]]) if tuple(point) in boundary_set]
        #bnd_mesh = nodes[bnd]
        nodes_dict = {tuple(node[:2]): node[2] for node in nodes}
        bnd_mesh = np.array([list(point) + [nodes_dict[tuple(point)]] for point in boundary_points])
    
    
    if Slab.Stress_drop_logic:
        fact_mu_z=np.zeros((Slab.N_scaling,len(barycenters_all[:,1])))
        for j in range(Slab.N_scaling):
            V1 = -(Slab.gamma2[0, j] + 2 * Slab.gamma1[0, j]) / (Slab.gamma1[0, j] + Slab.gamma2[0, j])
            V2 = -(Slab.gamma1[0, j] + 2 * Slab.gamma2[0, j]) / (Slab.gamma1[0, j] + Slab.gamma2[0, j])
            exponent = (Slab.gamma1[0, j] + Slab.gamma2[0, j]) / (Slab.gamma1[0, j] * V1 + Slab.gamma2[0, j] * V2 - Slab.gamma1[0, j] - Slab.gamma2[0, j])
    
            if not Slab.Rigidity_file_logic:
                mu_all, mu_BL, mu_bal = mesh.assign_rigidity(-1e-3 * barycenters_all[:, 2], Slab.Fact_rigidity, exponent)
                
            else:
                mu_all, mu_BL = mesh.assign_rigidity_from_file(-1e-3 * barycenters_all[:, 2], Slab.Rigidity_file)
                
            fact_mu_z[j, :] = mu_all / mu_BL
    else:
        if not Slab.Rigidity_file_logic:
            mu_all,_,_ = mesh.assign_rigidity(-1e-3 * barycenters_all[:, 2], Slab.Fact_rigidity)
        else:
            mu_all,_,_ = mesh.assign_rigidity_from_file(-1e-3 * barycenters_all[:, 2], Slab.Rigidity_file)
    
    #compute PDF for slip using coupling and rigidity
    K_all = mesh.coupling_pdf_CaA_function(1e-3 * barycenters_all[:, 2],-Slab.K_SL,-Slab.K_DL) #assign coupling
    SPDF_all = K_all / mu_all
    SPDF_all = SPDF_all / np.sum(SPDF_all)
    Slab.SPDF_all=SPDF_all
    
    
    name_filemu = os.path.join(main_dir,'config_files','Rigidity',f"mu_{Slab.zone_code}.dat")
    np.savetxt(name_filemu, mu_all, fmt="%.6f")
    
    
    #add other attributes to the Slab instance  
    Slab.bnd_mesh=bnd_mesh
    Slab.nodes=nodes
    Slab.cells=cells
    Slab.barycenters_all=barycenters_all
    Slab.hemisphere=hemisphere
    if Slab.Stress_drop_logic:   #FIXED by Antonio Scala
       Slab.fact_mu_z=fact_mu_z

def write_output(rupturing_areas, magnitudes, index_magnitude, spdf_all, name_scaling, area_cells, namefolder):
    print('Writing Output')
    
    for i , index in enumerate(index_magnitude):
        print(f'Magnitude bin # {index} - Mw={magnitudes[index]:.4f}')
        Mo = 10 ** (1.5 * magnitudes[index] + 9.1)
        folder_magnitude = f"{magnitudes[index]:6.4f}".replace('.', '_', 1)
        
        os.makedirs(folder_magnitude, exist_ok=True)
        os.chdir(folder_magnitude)
        
        for j, name_sc in enumerate(name_scaling):
            os.makedirs(name_sc, exist_ok=True)
            os.chdir(name_sc)
            
            for l, event in enumerate(rupturing_areas[i][j]):
                if event["true"]:
                    index_baryc = event["cell"][0]+1
                    
                    # QuakeArea file
                    filename = f'QuakeArea_{index_baryc:05d}.dat'
                    with open(filename, 'w') as fid:
                        fid.write('\n'.join(map(str, np.array(event["cell"])+1)) + '\n')
                    
                    # Slip_PDF file
                    spdf = spdf_all[event["cell"]]
                    spdf = spdf / sum(spdf)
                    filename = f'Slip_PDF_{index_baryc:05d}.dat'
                    with open(filename, 'w') as fid:
                        fid.write('\n'.join(f'{x:.6f}' for x in spdf) + '\n')
                    
                    # mu_Slip_aux file
                    area_event = sum(area_cells[event["cell"]]) * 1e6
                    slip2file = Mo / area_event
                    filename = f'mu_Slip_aux_{index_baryc:05d}.dat'
                    with open(filename, 'w') as fid:
                        fid.write(f'{slip2file:.4f}\n')
            
            os.chdir('..')
        os.chdir('..')
        shutil.move(folder_magnitude, os.path.join('..', namefolder, folder_magnitude))
        #os.rename(folder_magnitude, os.path.join('..', namefolder,folder_magnitude))



def generate_foldertree_slip(namefolder_slip, index_magnitude, Magnitude, Name_scaling):
    # Change to the parent directory
    #os.chdir('..')
    
    # Change to the specified folder
    #os.chdir(namefolder_slip)
    
    # Create main directories

    os.makedirs(os.path.join(main_dir,namefolder_slip,'homogeneous_mu'), exist_ok=True)
    os.makedirs(os.path.join(main_dir,namefolder_slip,'variable_mu'), exist_ok=True)
    
    for i in index_magnitude:
        folder_magnitude = f"{Magnitude[i]:6.4f}".replace('.', '_', 1)
        
        # Create subdirectories for each magnitude
        os.makedirs(os.path.join(main_dir,namefolder_slip,'homogeneous_mu',folder_magnitude), exist_ok=True)
        os.makedirs(os.path.join(main_dir,namefolder_slip,'variable_mu',folder_magnitude), exist_ok=True)
        
        for name in Name_scaling:
            os.makedirs(os.path.join(main_dir,namefolder_slip,'homogeneous_mu',folder_magnitude,name), exist_ok=True)
            os.makedirs(os.path.join(main_dir,namefolder_slip,'variable_mu',folder_magnitude,name), exist_ok=True)

    # Change back to the previous directory
    #os.chdir('..')






def read_file_lines(filename):
    with open(filename, 'r') as file:
        return file.read().splitlines()

def create_directory(path):
    if not os.path.exists(path):
        os.makedirs(path)

def run_command(command):
    subprocess.run(command, shell=True)
    
def run_command_Win(command,output_file):
    with open(output_file,'w') as outfile:
        subprocess.run([command, "input=param.dat"],stdout=outfile,text=True)

def copy_file(src, dst):
    if os.path.isfile(src):
        shutil.copy(src, dst)
        
def run_homo(slab):
    # Read input parameters from the slab instance

    event = slab.namefolder
    event_out = slab.namefolder_slip
    zone_code = slab.zone_code
    numb_stoch = slab.numb_stoch
    variable_mu = slab.variable_mu
    rigidity = 'homogeneous_mu'


    # Generate a list of magnitudes
    magnitudes = slab.get_magnitudes()
    classes_scaling=slab.Name_scaling
    folder_in = main_dir #os.getcwd()

    # Main processing loops
    for j, mw in enumerate(magnitudes):
        
        mw_string = f"{mw:6.4f}".replace('.', '_', 1)

        # Loop over each scaling class
        for i, cl in enumerate(classes_scaling):
            
            folder_out =os.path.join(main_dir,event_out,rigidity,mw_string,cl)# f"{folder_in}/{event_out}/{rigidity}/{mw_string}/{cl}"
            
            folder = os.path.join(main_dir,event,mw_string,cl) #f"{folder_in}/{event}/{mw_string}/{cl}"
            
            matrix_string = os.path.join(main_dir,'config_files','Matrix_distances',f"{zone_code}_matrix_distance.bin")#f"{folder_in}/config_files/Matrix_distances/{zone_code}_matrix_distance.bin\n"
            
            input_magnitude = f"magnitude={mw}\n"
            with open(os.path.join(folder_out,'input_magnitude'), 'w') as file:
                file.write(input_magnitude)

            
            with open(os.path.join(folder_out,'matrix_string.txt'), 'w') as file:
                file.write(matrix_string)

            #print(folder_out)

            # Copy necessary files
            shutil.copy(os.path.join(main_dir,'bin','k223d.x'), folder_out)
            #shutil.copy(os.path.join(main_dir,'input_magnitude'), folder_out)
            shutil.copy(os.path.join(main_dir,'param_zone.dat'), folder_out)
            #shutil.copy(os.path.join(main_dir,'matrix_string.txt'), folder_out)
            shutil.copy(os.path.join(main_dir,'config_files','Mesh',f'{zone_code}_mesh_15km.inp'), folder_out)
            shutil.copy(os.path.join(main_dir,'config_files','Parameters','param_homo.dat'), os.path.join(folder_out,'param.dat'))
        
            for quake_area_file in os.listdir(folder):

                if quake_area_file.startswith('QuakeArea') or quake_area_file.startswith('mu_Slip_aux'):
                    shutil.copy(os.path.join(folder, quake_area_file), folder_out)

            quake_area_files = [f for f in os.listdir(folder_out) if f.startswith('QuakeArea')]
            num_scenario = len(quake_area_files)
            with open(os.path.join(folder_out,'index_file.dat'), 'a') as file:
                file.write(f"{num_scenario} {numb_stoch}\n")

            for l in range(num_scenario):
                file_event = quake_area_files[l]
                eventid = file_event[-9:-4]

                for i in range(1, numb_stoch + 1):
                    if i < 10:
                        string4file = f"00{i}"
                    elif i < 100:
                        string4file = f"0{i}"
                    else:
                        string4file = str(i)

                    if i <= 2:
                        numb_gauss = 1
                    elif 2 < i <= 4:
                        numb_gauss = 2
                    else:
                        numb_gauss = 3

                    indexd = f"{eventid}_{string4file}"
                    with open(os.path.join(folder_out,'index_file.dat'), 'a') as file:
                        file.write(f"{indexd} {numb_gauss}\n")

            #shutil.move('index_file.dat', folder_out)
            os.chdir(folder_out)
            path_ex=os.path.join(folder_out,'k223d.x')
            current_os = platform.system()
            if current_os== 'Windows':
                print(folder_out)
                run_command_Win('k223d.x','output_file.log')
            else:
                print(folder_out)
                run_command('./k223d.x input=param.dat > output_file.log')
            
            #run_command(f'{path_ex} input=param.dat > output_file.txt')
            for file in quake_area_files:
                os.remove(os.path.join(folder_out,file))
            for file in os.listdir(folder_out):
                if file.startswith(('mu','Slip_PDF','Seed','in','param')):
                    os.remove(os.path.join(folder_out,file))
                if file.endswith(('txt','out','vtk','bin','inp','x')):
                    os.remove(os.path.join(folder_out,file))
            os.chdir(main_dir)
    
    if variable_mu == 0:
        create_directory(os.path.join(main_dir,'input'))
        create_directory(os.path.join(main_dir,'output'))
        shutil.move(os.path.join(main_dir,event_out), os.path.join(main_dir,'output'))
        shutil.move(os.path.join(main_dir,event), os.path.join(main_dir,'input'))
        #for file in [ 'input_magnitude']:
            #os.remove(os.path.join(main_dir,file))

        for file in os.listdir(main_dir):
            if file.endswith('.txt') or file.endswith('.dat'):
                os.remove(os.path.join(main_dir,file))
    

def run_var(slab):

    # Read input parameters from the file

    event = slab.namefolder
    event_out = slab.namefolder_slip
    zone_code = slab.zone_code
    numb_stoch = slab.numb_stoch
    variable_mu = slab.variable_mu
    rigidity = 'variable_mu'

    # Generate a list of magnitudes
    magnitudes = slab.get_magnitudes()
    classes_scaling=slab.Name_scaling
    
    folder_in = os.getcwd()

    # Main processing loops
    for j, mw in enumerate(magnitudes):
        
        mw_string = f"{mw:6.4f}".replace('.', '_', 1)

        # Loop over each scaling class
        for i, cl in enumerate(classes_scaling):

            folder = os.path.join(main_dir,event,mw_string,cl)#f"{folder_in}/{event}/{mw_string}/{cl}"
            folder_seed = os.path.join(main_dir,event_out,'homogeneous_mu',mw_string,cl)#f"{folder_in}/{event_out}/homogeneous_mu/{mw_string}/{cl}"
            folder_out = os.path.join(main_dir,event_out,rigidity,mw_string,cl)# f"{folder_in}/{event_out}/{rigidity}/{mw_string}/{cl}"

            os.makedirs(folder_out, exist_ok=True)

            with open(os.path.join(folder_out,'input_magnitude'), 'w') as f:
                f.write(f"magnitude={mw}\n")

            matrix_string = os.path.join(main_dir,'config_files','Matrix_distances',f"{zone_code}_matrix_distance.bin")#f"{main_dir}/config_files/Matrix_distances/{zone_code}_matrix_distance.bin"
            with open(os.path.join(folder_out,'matrix_string.txt'), 'w') as f:
                f.write(matrix_string + "\n")


             # Copy necessary files
            shutil.copy(os.path.join(main_dir,'bin','k223d.x'), folder_out)
            #shutil.copy(os.path.join(main_dir,'input_magnitude'), folder_out)
            shutil.copy(os.path.join(main_dir,'param_zone.dat'), folder_out)
            #shutil.copy(os.path.join(main_dir,'matrix_string.txt'), folder_out)
            shutil.copy(os.path.join(main_dir,'config_files','Mesh',f'{zone_code}_mesh_15km.inp'), folder_out)
            shutil.copy(os.path.join(main_dir,'config_files','Parameters','param_var.dat'), os.path.join(folder_out,'param.dat'))
            shutil.copy(os.path.join(main_dir,'config_files','Rigidity',f"mu_{zone_code}.dat"), folder_out)

            for quake_area_file in os.listdir(folder):

                if quake_area_file.startswith('QuakeArea') or quake_area_file.startswith('mu_Slip_aux') or quake_area_file.startswith('Slip_PDF'):
                    shutil.copy(os.path.join(folder, quake_area_file), folder_out)

       
            quake_area_files = [f for f in os.listdir(folder_out) if f.startswith('QuakeArea')]
            num_scenario = len(quake_area_files)

            with open(os.path.join(folder_out,'index_file.dat'), 'a') as index_file:
                index_file.write(f"{num_scenario} {numb_stoch}\n")

            for l in range(num_scenario):
                file_event = quake_area_files[l]
                eventid = file_event[-9:-4]

                for i in range(1, numb_stoch + 1):
                        if i < 10:
                            string4file = f"00{i}"
                        elif i < 100:
                            string4file = f"0{i}"
                        else:
                            string4file = str(i)

                        if i <= 2:
                            numb_gauss = 1
                        elif 2 < i <= 4:
                            numb_gauss = 2
                        else:
                            numb_gauss = 3

                        indexd = f"{eventid}_{string4file}"
                        with open(os.path.join(folder_out,'index_file.dat'), 'a') as file:
                            file.write(f"{indexd} {numb_gauss}\n")

            #shutil.move('index_file.dat', folder_out)
            os.chdir(folder_out)
            current_os = platform.system()
            if current_os== 'Windows':
                print(folder_out)
                run_command_Win('k223d.x','output_file.log')
            else:
                print(folder_out)
                run_command('./k223d.x input=param.dat > output_file.txt')

            for file in quake_area_files:
                os.remove(os.path.join(folder_out,file))
            for file in os.listdir(folder_out):
                if file.startswith(('mu','Slip_PDF','Seed','in','param')):
                    os.remove(os.path.join(folder_out,file))
                if file.endswith(('txt','out','vtk','bin','inp','x')):
                    os.remove(os.path.join(folder_out,file))
            os.chdir(main_dir)



    if variable_mu == 1:
        create_directory(os.path.join(main_dir,'input'))
        create_directory(os.path.join(main_dir,'output'))
        shutil.move(os.path.join(main_dir,event_out), os.path.join(main_dir,'output'))
        shutil.move(os.path.join(main_dir,event), os.path.join(main_dir,'input'))
        
        #for file in [ 'input_magnitude']:
            #os.remove(os.path.join(main_dir,file))
        for file in os.listdir(main_dir):
            if file.endswith('.txt') or file.endswith('.dat'):
                os.remove(os.path.join(main_dir,file))

