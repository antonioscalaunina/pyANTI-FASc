"""
June 2024

File to handle mesh files and compute mesh and mesh' elements attributes
"""

#Import libraries to be used
import json
import os
import time
import numpy as np
from tqdm import tqdm
from datetime import datetime
import sys

import scipy
from scipy.linalg import lstsq
from scipy.spatial import Delaunay
from scipy.spatial import KDTree

import utm
from pyproj import Transformer
from geopy.distance import geodesic
import geopy.distance as distance




###----------------------------------------------------------###
def read_nodes_cells(slab_name):
    
    print("Great! You already have nodes and cells, I'm just writing\n")
    
    nodes_filename = f"../utils/sz_slabs/{slab_name}/subfaults/{slab_name}_mesh_nodes.dat"
    if os.path.isfile(nodes_filename):
        nodes = np.loadtxt(nodes_filename)
    else:
        print('ERROR: Nodes file not in the sz database. Please check the name of slab or ensure your data is in the correct folder')
        return None, None
    
    # Adjusting longitude values
    nodes[nodes[:, 1] > 180, 1] -= 360
    
    faces_filename = f"../utils/sz_slabs/{slab_name}/subfaults/{slab_name}_mesh_faces.dat"
    if os.path.isfile(faces_filename):
        cells = np.loadtxt(faces_filename)
    else:
        print('ERROR: Faces file not in the sz database. Please check the name of slab or ensure your data is in the correct folder')
        return None, None
    
    return nodes, cells


###----------------------------------------------------------###

def fprintf(fid, format_str, *args):

    fid.write(format_str % args)

###----------------------------------------------------------###
# Function to extract nodes and create a unique list of nodes
def extract_nodes_and_cells(geojson_data):
    
    print("Great! I really love to create mesh files from GeoJSON format\n")
    
    nodes = []
    cells = []
    nodes_dict = {}  # Dictionary to store unique nodes (id_vertex)
    
    # Iterate over each feature in the geojson file
    for feature in geojson_data['features']:
        # Extract the vertices of the triangle and convert IDs to integers
        vertex1_id = int(feature['properties']['idvertex1'])  # Convert to int to remove leading zeros
        vertex1 = (vertex1_id, feature['properties']['lat1'], feature['properties']['lon1'], feature['properties']['depth1'])
        
        vertex2_id = int(feature['properties']['idvertex2'])  # Convert to int to remove leading zeros
        vertex2 = (vertex2_id, feature['properties']['lat2'], feature['properties']['lon2'], feature['properties']['depth2'])
        
        vertex3_id = int(feature['properties']['idvertex3'])  # Convert to int to remove leading zeros
        vertex3 = (vertex3_id, feature['properties']['lat3'], feature['properties']['lon3'], feature['properties']['depth3'])

        # Add the vertices to the dictionary if they are not already present
        for vertex in [vertex1, vertex2, vertex3]:
            id_vertex = vertex[0]
            if id_vertex not in nodes_dict:
                # Use the vertex ID directly from the JSON file
                nodes_dict[id_vertex] = {
                    'idx': id_vertex,  # Store the original ID
                    'lat': vertex[1],
                    'lon': vertex[2],
                    'depth': vertex[3]
                }
        
        # Create the cell with the original vertex IDs from the JSON file
        cell = [len(cells) + 1, vertex1_id, vertex2_id, vertex3_id]  # Cell number and vertex IDs
        cells.append(cell)

    # Convert the nodes dictionary to a list
    nodes = list(nodes_dict.values())
    
    # Sort nodes by their original ID (numerically) to ensure correct output order
    nodes.sort(key=lambda n: n['idx'])

    # Create a mapping from original ID to new index
    id_to_index = {node['idx']: idx + 1 for idx, node in enumerate(nodes)}

    # Update cells to use the new indices
    for cell in cells:
        cell[1] = id_to_index[cell[1]]  # Update vertex1 index
        cell[2] = id_to_index[cell[2]]  # Update vertex2 index
        cell[3] = id_to_index[cell[3]]  # Update vertex3 index

    return nodes, cells


def geojson2mesh(config_file):
    
    with open(config_file) as fid:
        Param = json.load(fid)
    
    #slab_acronym = Param['acronym'] #Slab
    nome_faglia = Param['zone_name'] #Slab
    
    # Step 1: Get the name of the fault (faglia) from user input
    # nome_faglia = input("Insert the name of the fault: ").strip()

    # Step 2: Define the directory structure
    base_dir = '../utils/sz_slabs'
    fault_dir = os.path.join(base_dir, nome_faglia)
    subfault_dir = os.path.join(fault_dir, 'subfaults')

    # Create directories if they don't exist
    os.makedirs(subfault_dir, exist_ok=True)

    # Step 3: Load the GeoJSON file (use the input name to locate it)
    geojson_filename = os.path.join(base_dir,f'{nome_faglia}_mesh.json')
    print(geojson_filename)
    if os.path.isfile(geojson_filename):
        with open(geojson_filename) as f:
            geojson_data = json.load(f)
    else:
        print('ERROR: Mesh in GeoJSON format does not exist! Please check option in input.json and zone/file names')
        sys.exit(1)

    # Step 4: Extract nodes and cells
    nodes, cells = extract_nodes_and_cells(geojson_data)

    # Step 5: Save nodes and cells with the correct filenames inside the "subfaults" folder
    nodes_filename = os.path.join(subfault_dir, f'{nome_faglia}_mesh_nodes.dat')
    cells_filename = os.path.join(subfault_dir, f'{nome_faglia}_mesh_faces.dat')

    # Save nodes to a .dat file
    with open(nodes_filename, 'w') as f:
        # Comment the header if needed
        # f.write('# idx lat lon depth\n')
        for node in nodes:  # Nodes are sorted by idx
            f.write(f"{node['idx']} {node['lon']} {node['lat']} {node['depth']}\n")

    # Save cells to a .dat file
    with open(cells_filename, 'w') as f:
        # Comment the header if needed
        # f.write('# cell vertex1 vertex2 vertex3\n')
        for cell in cells:  # Cells are already numbered in order of appearance
            f.write(f"{cell[0]} {cell[1]} {cell[2]} {cell[3]}\n")

    print(f"Files {nodes_filename} and {cells_filename} created successfully inside subfaults!")

### ----------------------------------------------------------------------- ####

def faces_nodes_2_mesh_file(config_file):

    """
    Generate mesh input file from nodes and faces files

    """

    with open(config_file) as fid:
        Param = json.load(fid)
    
    slab_acronym = Param['acronym'] #Slab
    slab_name = Param['zone_name'] #Slab

    string1=scipy.io.loadmat('../utils/string1.mat')['string1'][0]
    string1[1][0] = f'cubit({os.getcwd()}):{datetime.now().strftime("%Y-%m-%d")}:'
    string2=scipy.io.loadmat('../utils/string2.mat')['string2'][0]
    string3=scipy.io.loadmat('../utils/string3.mat')['string3'][0]
    
    nodes, cells= read_nodes_cells(slab_name)
    namefile = f"../config_files/Mesh/{slab_acronym}_mesh_15km.inp"

    with open(namefile, 'w') as fid:
        for line in string1:
            fid.write(line[0] + '\n')

        for node in nodes:
            fprintf(fid, '%d, %20.10e, %20.10e, %15.6e\n', *node)

        for line in string2:
            fid.write(line[0] + '\n')

        for cell in cells:
            fprintf(fid, '%d, %8d, %8d, %8d\n', *cell)

        for line in string3:
            fid.write(line[0] + '\n')

    #compute distance mstrix
    #matrix_distance=matrix_distance_nolat2(nodes[:,1:4])

    #save matrix for further use
    #file_path = f'../config_files/Matrix_distances/{slab_acronym}_matrix_distance.bin'
    #with open(file_path, 'wb') as f:
        #f.write(matrix_distance.tobytes())


###----------------------------------------------------------###

def read_mesh_file(fid):
    """
    Read file mesh in the variable fid in inp format providing as nodes
    (Lon-Lat-Depth) and the cells (node indices of the vertices)
    """
    # Skip first 9 lines
    for _ in range(9):
        fid.readline()

    # Read nodes
    for i in range(100000):
        line = fid.readline()
        if line.startswith('*'):
            numb_nodes = i
            break

    # Skip next 2 lines
    for _ in range(2):
        fid.readline()

    # Read cells
    for i in range(100000):
        line = fid.readline()
        if line.startswith('*'):
            numb_cells = i
            break

    # Rewind file
    fid.seek(0)
    nodes = np.zeros((numb_nodes, 3))
    cells = np.zeros((numb_cells, 3), dtype=int)
    string1 = []
    string2 = []
    string3 = []

    # Read first 9 lines again into string1
    for _ in range(9):
        string1.append(fid.readline().strip())

    # Read nodes
    for i in range(numb_nodes):
        line = fid.readline().strip()
        temp = list(map(float, line.split(',')))
        nodes[i, :] = temp[1:4]

    # Read next 3 lines into string2
    for _ in range(3):
        string2.append(fid.readline().strip())

    # Read cells
    for i in range(numb_cells):
        line = fid.readline().strip()
        temp = list(map(int, line.split(',')))
        cells[i, :] = temp[1:4]

    # If there are more outputs, read remaining 44 lines into string3
    if len(fid.readlines()) > 0:
        fid.seek(0)
        for _ in range(numb_nodes + numb_cells + 14):
            fid.readline()  # Skip lines already read
        for _ in range(44):
            string3.append(fid.readline().strip())

    return nodes, cells, string1, string2, string3

###----------------------------------------------------------###

def compute_connectivity(EToV):
    """
    Compute connectivity matrix of the elements in EToV
    """
    Nfaces = 3
    K = EToV.shape[0]
    Nnodes = np.max(EToV)

    # Create list of all faces 1, then 2, & 3
    fnodes = np.vstack([EToV[:, [0, 1]], EToV[:, [1, 2]], EToV[:, [2, 0]]])
    fnodes = np.sort(fnodes, axis=1) - 1
    
    # Set up default element to element and Element to faces connectivity
    EToE = np.tile(np.arange(1, K + 1).reshape(K, 1), (1, Nfaces))
    EToF = np.tile(np.arange(1, Nfaces + 1), (K, 1))

    # Uniquely number each set of three faces by their node numbers
    id = fnodes[:, 0] * Nnodes + fnodes[:, 1] + 1
    spNodeToNode = np.column_stack([id, np.arange(1, Nfaces * K + 1), EToE.T.ravel(), EToF.T.ravel()])
    
    EToE = np.zeros((K, Nfaces))
    EToF = np.zeros((K, Nfaces))

    # Now we sort by global face number
    sorted_indices = np.argsort(spNodeToNode[:, 0])
    sorted_spNodeToNode = spNodeToNode[sorted_indices]
   
    # Find matches in the sorted face list
    match_indices = np.where(sorted_spNodeToNode[:-1, 0] == sorted_spNodeToNode[1:, 0])[0]
    
    # Make links reflexive
    matchL = np.vstack([sorted_spNodeToNode[match_indices], sorted_spNodeToNode[match_indices + 1]])
    matchR = np.vstack([sorted_spNodeToNode[match_indices + 1], sorted_spNodeToNode[match_indices]])
    
    # Insert matches
    EToE.T.flat[matchL[:, 1] - 1] = matchR[:, 2]
    EToF.T.flat[matchL[:, 1] - 1] = matchR[:, 3]

    return EToE, EToF

###----------------------------------------------------------###


def haversine_vectorized(lat1, lon1, lat2, lon2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees) in a vectorized manner.
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])

    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a)) 
    km = 6371 * c # Radius of earth in kilometers. Use 3956 for miles. Determines return value units.
    return km * 1000 # return in meters

def matrix_distance_nolat2(nodes,slab_acronym):
    """
    Compute the distance matrix for a set of nodes, considering both horizontal and vertical distances.

    Parameters:
        nodes (np.ndarray): A 2D array where each row represents a node with coordinates (lon, lat, depth).

    Returns:
        np.ndarray: A symmetric matrix of distances between the nodes.
    """
    n = len(nodes)
    matrix_distance = np.zeros((n, n), dtype=np.float32)

    # Extract coordinates and depths
    lons = nodes[:, 0]
    lats = nodes[:, 1]
    depths = nodes[:, 2]

    # Vectorized computation of haversine distances
    lon1, lon2 = np.meshgrid(lons, lons)
    lat1, lat2 = np.meshgrid(lats, lats)
    horizontal_distances = haversine_vectorized(lat1, lon1, lat2, lon2)

    # Vectorized computation of vertical distances
    depth1, depth2 = np.meshgrid(depths, depths)
    vertical_distances = np.abs(depth2 - depth1)

    # Combine horizontal and vertical distances
    matrix_distance = np.sqrt(horizontal_distances**2 + vertical_distances**2)

    print('Distance matrix computed!')

    #save matrix for further use
    file_path = f'../config_files/Matrix_distances/{slab_acronym}_matrix_distance.bin'
    with open(file_path, 'wb') as f:
        f.write(matrix_distance.tobytes())

    return matrix_distance

###-----------------------------------------------------------###

def matrix_distance_nolat(nodes):
    """
    Compute the distance matrix for a set of nodes, considering both horizontal and vertical distances.
    
    Parameters:
        nodes (np.ndarray): A 2D array where each row represents a node with coordinates (lon, lat, depth).
        
    Returns:
        np.ndarray: A symmetric matrix of distances between the nodes.
    """
    n = len(nodes)
    matrix_distance = np.zeros((n, n), dtype=np.float32)

    progress_bar = tqdm(total=n, desc="Computing distance matrix")
    
    for i in range(n):
        progress_bar.update(1)
        for j in range(i + 1, n):
            
            distance_wh = distance.distance((nodes[i,1], nodes[i,0]), (nodes[j,1], nodes[j,0])).meters
            matrix_distance[i, j] = np.sqrt(distance_wh**2 + (nodes[j, 2] - nodes[i, 2])**2)
            matrix_distance[j, i] = matrix_distance[i, j]
            
    progress_bar.close()
    print('Distance matrix computed!')
    
    return matrix_distance


###----------------------------------------------------------###

def compute_area(cells,matrix_distance):
    """
    Compute area for a set of elements (cells) and the total area

    Parameters:
        nodes (np.array):
        cells (np:array):
        matrix_distance (np:array); 
    Returns:
        
    
    """
    Area_cells=np.zeros(len(cells))
    #X,Y,UTMS,_=utm.from_latlon(self.nodes[:,1], self.nodes[:,0], force_zone_number=self.Merc_zone)
    
    for i, cell in enumerate(cells):
    
        #x,y,_,_=utm.from_latlon(nodes[cell-1,1], nodes[cell-1,0])
        #a=np.linalg.norm(np.array([x[0], y[0],nodes[cell[0]-1,2]] ) - np.array([x[1],y[1],nodes[cell[1]-1,2]] ))
        #b=np.linalg.norm(np.array([x[0], y[0],nodes[cell[0]-1,2]] ) - np.array([x[2],y[2],nodes[cell[2]-1,2]] ))
        #c=np.linalg.norm(np.array([x[2],y[2],nodes[cell[2]-1,2]] ) - np.array([x[1],y[1],nodes[cell[1]-1,2]] ))
        a=matrix_distance[cell[0]-1,cell[1]-1]
        b=matrix_distance[cell[0]-1,cell[2]-1]
        c=matrix_distance[cell[1]-1,cell[2]-1]

        s = (a + b + c) / 2
        Area_cells[i] = np.sqrt(s * (s - a) * (s - b) * (s - c))
    
    Area_tot=np.sum(Area_cells)
    return Area_tot, Area_cells




def find_barycenters(nodes, cells):
    """
    Calculate the barycenters of the given cells using their node coordinates.

    Parameters:
    nodes (numpy.ndarray): Array of node coordinates (Lon, Lat, Depth).
    cells (numpy.ndarray): Array of cells' nodes indices.

    Returns:
    numpy.ndarray: Array of barycenters for each cell.
    """

    # Adjust longitude values
    nodes[nodes[:, 0] < 0, 0] += 360

    # Initialize barycenters array
    Barycenters = np.zeros((cells.shape[0], 3))

    # Calculate barycenters
    for i in range(cells.shape[0]):
        Barycenters[i, 0] = np.sum(nodes[cells[i, :]-1, 0]) / 3
        Barycenters[i, 1] = np.sum(nodes[cells[i, :]-1, 1]) / 3
        Barycenters[i, 2] = np.sum(nodes[cells[i, :]-1, 2]) / 3

    # Adjust longitude values for barycenters
    Barycenters[Barycenters[:, 0] > 180, 0] -= 360

    return Barycenters



def assign_rigidity(depth, Fact_mu, exponent=None):
    """
    Assign rigidity variation
    """
    PREM_mu = np.array([
        [0, 26.5], [3, 26.5], [15, 43.9], [24.4, 67.8],
        [40, 67.6], [60, 67.3], [80, 64.6], [115, 64.2]
    ])

    a = 0.5631
    b = 0.0437
    mu = np.zeros(len(depth))
    mu_BL = np.zeros(len(depth))
    rigidity_PREM = np.zeros(len(depth))
    
    # GEIST & BILEK RIGIDITY
    # AVERAGE WITH PREM RIGIDITY
    for i in range(len(depth)):
        mu[i] = 10**(a + b * depth[i])
        mu_BL[i] = mu[i]
        for j in range(len(PREM_mu) - 1):
            if depth[i] >= PREM_mu[j, 0] and depth[i] <= PREM_mu[j + 1, 0]:
                rigidity_PREM[i] = linear_interp(depth[i], PREM_mu[j:j + 2, :])
                break
        
        if mu[i] >= rigidity_PREM[i]:
            if depth[i] > 10:
                mu[i] = rigidity_PREM[i]
                mu_BL[i] = rigidity_PREM[i]
        else:
            mu[i] = mu[i] + Fact_mu * (rigidity_PREM[i] - mu[i])
    
    mu_bal = None
    if exponent is not None:
        index = np.argmin(np.abs(mu_BL - rigidity_PREM))
        coeff = mu_BL[index]**(-1)
        mu_bal = coeff**exponent * mu_BL**(exponent + 1)
    
    return mu, mu_BL, mu_bal


def assign_rigidity_from_file(depth, namefile):
    '''
    Compute and assign depth-dependent rigidity to depths in the depth array from the file: namefile
    '''
    PREM_mu = np.array([
        [0, 26.5], [3, 26.5], [15, 43.9], [24.4, 67.8],
        [40, 67.6], [60, 67.3], [80, 64.6], [115, 64.2]
    ])

    a = 0.5631
    b = 0.0437

    mu_BL = np.zeros(len(depth))
    rigidity_PREM = np.zeros(len(depth))
    
    # GEIST & BILEK RIGIDITY
    # AVERAGE WITH PREM RIGIDITY
    for i in range(len(depth)):
        mu_BL[i] = 10**(a + b * depth[i])
        for j in range(len(PREM_mu) - 1):
            if depth[i] >= PREM_mu[j, 0] and depth[i] <= PREM_mu[j + 1, 0]:
                rigidity_PREM[i] = linear_interp(depth[i], PREM_mu[j:j + 2, :])
                break
        
        if mu_BL[i] >= rigidity_PREM[i]:
            if depth[i] > 10:
                mu_BL[i] = rigidity_PREM[i]
    
    if not os.path.exists(f'../config_files/Parameters/{namefile}'):
        print('Error: File of rigidity variation not available: PLEASE CHECK!')
        return None, None
    else:
        table_mu = np.genfromtxt(f'../config_files/Parameters/{namefile}', delimiter=',', skip_header=1)
        
        mu = np.zeros(len(depth))
        for i in range(len(depth)):
            for j in range(len(table_mu) - 1):
                if depth[i] >= table_mu[j, 0] and depth[i] <= table_mu[j + 1, 0]:
                    mu[i] = linear_interp(depth[i], table_mu[j:j + 2, :])
                    break

    mu[depth < table_mu[0, 0]] = table_mu[0, 1]
    mu[depth > table_mu[-1, 0]] = table_mu[-1, 1]

    return mu, mu_BL




def linear_interp(depth, PREM):
    m = (np.log10(PREM[1, 1]) - np.log10(PREM[0, 1])) / (PREM[1, 0] - PREM[0, 0])
    q = np.log10(PREM[0, 1]) - m * PREM[0, 0]
    mu = m * depth + q
    return 10**mu




def coupling_pdf_CaA_function(depth, up_l, low_l):
    '''
    Compute and assign coupling factor K to depths in array depth
    '''
    
    
    V1 = [low_l, 1]
    P1 = [low_l-10, 0]

    # Calculate coefficients for the parabolic section
    B1 = np.array([
        [-V1[0]**2, 1, V1[1]],
        [P1[0]**2 - P1[0] * 2 * V1[0], 1, P1[1]]
    ])
    B1_aux = lstsq(B1[:, :-1], B1[:, -1])[0]

    a = B1_aux[0]
    b = -a * 2 * V1[0]
    c = B1_aux[1]

    # Calculate coefficients for the 7th power polynomial section
    A_pow = np.array([
        [(-13)**7, 1, 1],
        [(-3.9816)**7, 1, 0]
    ])
    coeff_pow = lstsq(A_pow[:, :-1], A_pow[:, -1])[0]

    a_pow = coeff_pow[0]
    b_pow = coeff_pow[1]

    # Calculate coefficients for the shallow parabolic link
    V2 = [up_l, 1]
    P2 = [up_l+2.5, a_pow * (up_l+2.5)**7 + b_pow]

    B2 = np.array([
        [-V2[0]**2, 1, V2[1]],
        [P2[0]**2 - P2[0] * 2 * V2[0], 1, P2[1]]
    ])
    B2_aux = lstsq(B2[:, :-1], B2[:, -1])[0]

    a_d = B2_aux[0]
    b_d = -a_d * 2 * V2[0]
    c_d = B2_aux[1]


    
    # Calculate K values based on depth
    K = np.zeros(len(depth))
    for i in range(len(depth)):
        if depth[i] <= low_l:
            K[i] = a * depth[i]**2 + b * depth[i] + c
        elif low_l < depth[i] <= up_l:
            K[i] = 1
        elif up_l < depth[i] < up_l+2.5:
            K[i] = a_d * depth[i]**2 + b_d * depth[i] + c_d
        else:
            K[i] = a_pow * depth[i]**7 + b_pow

        if K[i] < 0:
            K[i] = 0

    return K


####################################################################

def get_hemisphere(lat):
    """
    Set hemisphere for lat,lon to x,y conversion based on mesh' latitudes
    """
    lat_mean=np.mean(lat)
    if lat_mean>0:
        hemisphere='N'
    else:
        hemisphere='S'
    
    return hemisphere


def get_utm_epsg(zone, hemisphere='N'):
    if hemisphere == 'N':
        epsg_code = 32600 + zone
    elif hemisphere == 'S':
        epsg_code = 32700 + zone
    else:
        raise ValueError("Hemisphere must be 'N' or 'S'")
    return epsg_code


def ll2xy(lon,lat,zone,hemisphere):
    
    source_crs = "EPSG:4326"
    target_crs = get_utm_epsg(zone, hemisphere)
    # Create a Transformer object
    transformer = Transformer.from_crs(source_crs, target_crs, always_xy=True)
    x, y = transformer.transform(lon, lat)
    return x,y



def from_latlon(lat, lon, zone):
    # Check if inputs are single values or arrays
    if np.isscalar(lat) and np.isscalar(lon):
        # Single value input
        lon = np.array([lon])
        lat = np.array([lat])
        single_value = True

        lon[lon > 180] -= 360

        x, y, _, _ = utm.from_latlon(lat, lon, force_zone_number=zone)
        
        if lat[0]<0:
            y-=10000000 # Convert to Northern Hemisphere coordinates

        return x[0], y[0]
    
    else:
        # Array input
        lon = np.asarray(lon)
        lat = np.asarray(lat)
        single_value = False

        # Adjust longitude values
        lon[lon > 180] -= 360

        # Initialize northing and easting arrays
        x = np.zeros(len(lon))
        y = np.zeros(len(lat))

        # Positive sign latitudes
        if len(lat[lat > 0])>0:
            x[lat > 0], y[lat > 0], _, _ = utm.from_latlon(lat[lat > 0], lon[lat > 0], force_zone_number=zone)

        # Negative sign latitudes
        if len(lat[lat > 0])<0:
            x[lat < 0], y1, _, _ = utm.from_latlon(lat[lat < 0], lon[lat < 0], force_zone_number=zone)
            y[lat < 0] = y1 - 10000000  # Convert to Northern Hemisphere coordinates


        return x, y
############################################################################################
