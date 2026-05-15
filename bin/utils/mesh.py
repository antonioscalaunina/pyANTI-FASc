"""
June 2024

File to handle mesh files and compute mesh and mesh' elements attributes
"""

#Import libraries to be used
import json
import csv
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
from pyproj import CRS, Transformer
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

    coords_to_id = {}   # (lat, lon, depth) -> node_id
    nodes_dict = {}     # node_id -> node
    next_node_id = 1

    def get_vertex(feature, lat_key, lon_key, depth_key):
        nonlocal next_node_id

        lat = feature['properties'][lat_key]
        lon = feature['properties'][lon_key]
        depth = feature['properties'][depth_key]

        coords = (
            round(lat, 8),
            round(lon, 8),
            round(depth, 6)
        )

        if coords in coords_to_id:
            return coords_to_id[coords]

        node_id = next_node_id
        next_node_id += 1

        coords_to_id[coords] = node_id

        nodes_dict[node_id] = {
            'idx': node_id,
            'lat': lat,
            'lon': lon,
            'depth': depth
        }

        return node_id

    for feature in geojson_data['features']:

        v1 = get_vertex(feature, 'lat1', 'lon1', 'depth1')
        v2 = get_vertex(feature, 'lat2', 'lon2', 'depth2')
        v3 = get_vertex(feature, 'lat3', 'lon3', 'depth3')

        # evita triangoli degeneri
        if len({v1, v2, v3}) < 3:
            print("⚠️ Degenerate triangle skipped")
            continue

        cell = [len(cells) + 1, v1, v2, v3]
        cells.append(cell)

    nodes = list(nodes_dict.values())
    nodes.sort(key=lambda n: n['idx'])

    return nodes, cells

###------------------------------------------------------------###
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
    geojson_filename = os.path.join(base_dir,f'{nome_faglia}_mesh.geojson')
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
    nodes, cells= read_nodes_cells(slab_name)
    n_cells = len(cells)

    # Modification to allow user to define rake either by a single value or through a file (2026/04/01)
    # If either value is not valid, or file is not found, or no rake is defined standard Rake=90° is used 
    
    if "rake" in Param:
        rake_file = Param["rake"]
        if isinstance(rake_file,str):
            try:
                # Open the CSV file
                with open(rake_file, newline='', encoding="utf-8") as f:
                    reader = csv.reader(f)
                    original_headers = next(reader)

                    # Normalize headers: lowercase and strip spaces
                    headers = [h.strip().lower() for h in original_headers]

                    # Create DictReader with normalized headers
                    f.seek(0)
                    dict_reader = csv.DictReader(f, fieldnames=headers)
                    next(dict_reader)  # skip the original header row

                    # Store all values from the 'rake' column into Slab.rakes
                    rakes = [row['rake'] for row in dict_reader]

            except FileNotFoundError:
                    # Warn if the file is missing
                    print(
                        f"Rake file '{rake_file}' not found: check path and filename. A standard rake=90° will be used",
                    )
        elif isinstance(rake_file,(int,float)):
            if  -180 <= rake_file <= 180:
                rakes = [float(rake_file)] * n_cells
            else:
                print(
                    f"Invalid rake value ({rake_source}). Must be between -180 and 180. "
                    "A standard rake=90° will be used"
                )

        else:
            print(
                f"Unsupported rake input type ({type(rake_file)}). "
                "A standard rake=90° will be used"
            )
    else:
        print("Rake not defined: a standard rake=90° will be used")

# END OF MODIFICATION FOR RAKE!


    string1=scipy.io.loadmat('../utils/string1.mat')['string1'][0]
    string1[1][0] = f'cubit({os.getcwd()}):{datetime.now().strftime("%Y-%m-%d")}:'
    string2=scipy.io.loadmat('../utils/string2.mat')['string2'][0]
    string3=scipy.io.loadmat('../utils/string3.mat')['string3'][0]
    


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
        
        #If a rake is defined the mesh will contain this info (2026/04/01)

        if 'rakes' in locals():
            if len(rakes) == len(cells):
                # 1. write first two lines of string3
                for line in string3[:2]:
                    fid.write(line[0] + '\n')

                # 2. write RAKE field header
                fid.write("** FIELD VARIABLE 1 = RAKE (degrees)\n")
                fid.write("*INITIAL CONDITIONS, TYPE=FIELD, VARIABLE=1\n")

                # 3. write each rake value with its index (1-based)
                for idx, rake_val in enumerate(rakes, start=1):
                    fid.write(f"{idx}, {rake_val}\n")

                # 4. write the remaining lines of string3
                for line in string3[2:]:
                    fid.write(line[0] + '\n')

            else:
                # rake exists but dimensions mismatch
                print("Warning: RAKE and cells dimensions do not match. Using standard rake=90°. Writing string3")
                for line in string3:
                    fid.write(line[0] + '\n')

        else:
            # rake not defined at all
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
    #print(Area_tot,Area_cells)
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
    
    rigidity_file = namefile

    if rigidity_file is None:
        print("Error: Rigidity_file is not defined in input.json")
        return None, None
    
    # Accept:
    # 1. absolute paths
    # 2. relative paths already valid from current working directory
    # 3. simple filenames to be searched in ../config_files/Rigidity/
    if not os.path.isabs(rigidity_file):
        if os.path.isfile(rigidity_file):
            pass
        else:
            rigidity_file = os.path.join("..", "config_files", "Rigidity", rigidity_file)
    
    if not os.path.exists(rigidity_file):
        print(f"Error: File of rigidity variation not available: {namefile}")
        return None, None
    
    table_mu = np.genfromtxt(
        rigidity_file,
        delimiter=",",
        skip_header=1
    )
    
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

# ==========================================================
# Coordinate transformations
# ==========================================================

def local_projection(lon0, lat0):
    """
    Create local azimuthal equidistant projection centered on (lon0, lat0)
    """
    crs_geo = CRS.from_epsg(4326)
    crs_local = CRS.from_proj4(
        f"+proj=aeqd +lat_0={lat0} +lon_0={lon0} +units=m +datum=WGS84"
    )

    fwd = Transformer.from_crs(crs_geo, crs_local, always_xy=True)
    inv = Transformer.from_crs(crs_local, crs_geo, always_xy=True)

    return fwd, inv


# ==========================================================
# Fault plane mesh generation (metric space)
# ==========================================================

def fault_plane_mesh(
    center_xyz,
    length,
    width,
    strike_deg,
    dip_deg,
    target_area
):
    """
    Build a triangular mesh of a rectangular fault plane in metric coordinates.
    """

    strike = np.deg2rad(strike_deg)
    dip = np.deg2rad(dip_deg)

    # Unit vector along strike (horizontal)
    strike_vec = np.array([
        np.sin(strike),
        np.cos(strike),
        0.0
    ])

    # Unit vector along dip (downward)
    dip_vec = np.array([
        np.cos(strike) * np.cos(dip),
       -np.sin(strike) * np.cos(dip),
       -np.sin(dip)
    ])

    # Grid resolution
    cell_size = np.sqrt(2 * target_area)
    n_strike = max(1, int(np.ceil(length / cell_size)))
    n_dip = max(1, int(np.ceil(width / cell_size)))

    s = np.linspace(-length / 2, length / 2, n_strike + 1)
    d = np.linspace(-width / 2, width / 2, n_dip + 1)

    # Vertices
    vertices = []
    for di in d:
        for si in s:
            p = center_xyz + si * strike_vec + di * dip_vec
            vertices.append(p)

    vertices = np.array(vertices)

    # Triangles
    triangles = []

    def idx(i, j):
        return i * (n_strike + 1) + j

    for i in range(n_dip):
        for j in range(n_strike):
            v0 = idx(i, j)
            v1 = idx(i, j + 1)
            v2 = idx(i + 1, j)
            v3 = idx(i + 1, j + 1)

            triangles.append([v0, v1, v2])
            triangles.append([v1, v3, v2])

    return vertices, np.array(triangles)


# ==========================================================
# Create triangles in lonlat (geographic coordinates)
# ==========================================================

def triangles_to_lonlat(vertices, triangles, transformer_inv):
    """
    Convert mesh triangles from local coordinates to geographic coordinates.
    """

    triangles_lonlat = []
    triangles_depth = []

    for tri in triangles:
        tri_ll = []
        tri_z = []

        for idx in tri:
            x, y, z = vertices[idx]
            lon, lat = transformer_inv.transform(x, y)

            tri_ll.append((lon, lat))
            tri_z.append(z)  # depth negative downward

        triangles_lonlat.append(tri_ll)
        triangles_depth.append(tri_z)

    return triangles_lonlat, triangles_depth

# ==========================================================
# Save extended geoJSON (QGIS format)
# ==========================================================

def save_extended_mesh_geojson(
    triangles_lonlat,
    triangles_depth,
    strike,
    dip,
    filename
):
    """
    Save mesh in an extended GeoJSON format similar to GRCF003_mesh.json.
    Geometry is 2D (lon, lat); depth is stored in properties.
    """

    features = []

    for i, (tri_ll, tri_z) in enumerate(
        zip(triangles_lonlat, triangles_depth), start=1
    ):
        (lon1, lat1), (lon2, lat2), (lon3, lat3) = tri_ll
        z1, z2, z3 = tri_z

        # Centroid
        lon_c = (lon1 + lon2 + lon3) / 3
        lat_c = (lat1 + lat2 + lat3) / 3
        depth_c = (z1 + z2 + z3) / 3

        feature = {
            "type": "Feature",
            "id": f"mesh.{i:06d}",
            "geometry": {
                "type": "Polygon",
                "coordinates": [[
                    [lon1, lat1],
                    [lon2, lat2],
                    [lon3, lat3],
                    [lon1, lat1]
                ]]
            },
            "properties": {
                "idtriangle": i,

                "lon1": lon1, "lat1": lat1, "depth1": z1,
                "lon2": lon2, "lat2": lat2, "depth2": z2,
                "lon3": lon3, "lat3": lat3, "depth3": z3,

                "lon_c": lon_c,
                "lat_c": lat_c,
                "depth_c": depth_c,

                "strike": strike,
                "dip": dip
            }
        }

        features.append(feature)

    geojson = {
        "type": "FeatureCollection",
        "features": features
    }

    with open(filename, "w") as f:
        json.dump(geojson, f, indent=2)

    print(f"Extended GeoJSON written to: {filename}")


def point2geojson(config_file):

    
    
    with open(config_file) as fid:
        Param = json.load(fid)
    zone_name=Param["zone_name"]
    filename = f"../utils/sz_slabs/{zone_name}_mesh.geojson"
    # ------------------------------
    # Input fault parameters
    if "lon_c" in Param:
        lon_val = Param["lon_c"]
        if -180.0 <= lon_val <= 180.0:
            lon0 = lon_val
        else:
            raise ValueError(f"Invalid longitude (lon_c = {lon_val}). Must be in the range [-180, 180]")
    else:
        raise KeyError(f"Longitude not defined check 'mesh_gen' or lon_c parameter")
    
    if "lat_c" in Param:
        lat_val = Param["lat_c"]
        if -90.0 <= lat_val <= 90.0:
            lat0 = lat_val
        else:
            raise ValueError(f"Invalid longitude (lat_c = {lat_val}). Must be in the range [-90, 90]")
    else:
        raise KeyError(f"Latitude not defined check 'mesh_gen' or lat_c parameter")
    

    if "depth_km" in Param:
        depth0 = Param['depth_km']  
        depth0 = 1000.0*depth0 # meters (positive downward)
    else:
        print("Depth not defined in the input.json file. A standard depth of 10 km is used")
        depth0 = 10000.0 #

    if "length_km" in Param:
        length = Param['length_km']  
        length = 1000.0*length # meters 
    else:
        print("Length not defined in the input.json file. A standard length of 50 km is used")
        length = 50000.0 #

    if "width_km" in Param:
        width = Param['width_km']  
        width = 1000.0*width # meters 
    else:
        print("Width not defined in the input.json file. A standard width of 25 km is used")
        width = 25000.0 #

    if "strike" in Param:
        strike_val = Param["strike"]
        if 0.0 <= strike_val <= 360.0:
            strike = strike_val
        else:
            raise ValueError(f"Invalid strike (strike = {strike_val}). Must be in the range [0, 360]")
    else:
        print("Strike not defined in the input.json file. A standard strike of 0° is used")
        strike=0.0

    if "dip" in Param:
        dip_val = Param["dip"]
        if 0.0 <= dip_val <= 90.0:
            dip = dip_val
        else:
            raise ValueError(f"Invalid dip (dip = {dip_val}). Must be in the range [0, 90]")
    else:
        print("Dip not defined in the input.json file. A standard dip of 30° is used")
        dip=30.0

    if "elem_size_km2" in Param:
        area = Param['elem_size_km2']  
        area = 1e6*area # meters^2 
    else:
        print("Average element size not defined in the input.json file. A standard area of 5 km^2 is used")
        area = 5.0e6 # m^2 (average triangle area)


    # ------------------------------
    # Build local projection
    fwd, inv = local_projection(lon0, lat0)

    x0, y0 = fwd.transform(lon0, lat0)
    center_xyz = np.array([x0, y0, -depth0])

    # ------------------------------
    # Generate mesh
    vertices, triangles = fault_plane_mesh(
        center_xyz,
        length,
        width,
        strike,
        dip,
        area
    )
    # ------------------------------
    # Check profondità (devono essere tutte <= 0)
    if np.any(vertices[:, 2] > 0):
        raise ValueError(
            "ERROR: some depths are positive (above the surface). "
            "Increase the central depth or adjust the fault parameters."
        )
    # ------------------------------
    # Plot
    #plot_fault_mesh(vertices, triangles)

    # ------------------------------
    # Export GeoJSON
    #save_mesh_geojson(vertices, triangles, inv, "fault_mesh.geojson")

    # -------------------------------------------------
    # Convert triangles to geographic coordinates
    triangles_lonlat, triangles_depth = triangles_to_lonlat(
        vertices,
        triangles,
        inv
    )

    # -------------------------------------------------
    # Save EXTENDED GeoJSON (GRCF-like)
    save_extended_mesh_geojson(
        triangles_lonlat,
        triangles_depth,
        strike,
        dip,
        filename
    )

