# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 15:28:14 2024

@author: ascal
"""

import json
import os

# Function to extract nodes and create a unique list of nodes
def extract_nodes_and_cells(geojson_data):
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

# Main code
if __name__ == "__main__":
    # Step 1: Get the name of the fault (faglia) from user input
    nome_faglia = input("Insert the name of the fault: ").strip()

    # Step 2: Define the directory structure
    base_dir = 'sz_slabs'
    fault_dir = os.path.join(base_dir, nome_faglia)
    subfault_dir = os.path.join(fault_dir, 'subfaults')

    # Create directories if they don't exist
    os.makedirs(subfault_dir, exist_ok=True)

    # Step 3: Load the GeoJSON file (use the input name to locate it)
    geojson_filename = f'{nome_faglia}_mesh.json'
    with open(geojson_filename) as f:
        geojson_data = json.load(f)

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


