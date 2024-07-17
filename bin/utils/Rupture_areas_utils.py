"""
June 2024

"""

#Import libraries to be used
import json
import os
import numpy as np
import scipy
from scipy.linalg import lstsq
from scipy.spatial import Delaunay
import utm

import alphashape
from shapely.geometry import Point, Polygon
from shapely.vectorized import contains
from geopy.distance import geodesic
import geopy.distance as distance

import utils.mesh as mesh


def find_closest_point_segment(segment, q):
    """
    Calculate the closest point of the segment to point(s) in q

    Parameters:
    segment (numpy.ndarray(2,2)): Array containing the coordinates (easting, northing) of the end-points of the fault segment
    q (numpy.ndarray(N,2)): Array of coordinates (easting, northing) of points of interest (e.g., barycenters, stations).

    Returns:
    p (numpy.ndarray (N,2)): Array containing the coordinates (easting, northing) of the closest point of the segment to each point in q
    dist (numpy.ndarray (N)): Array of distance between the closest point of the segment to each point in q
    """
    
    p1 = segment[0,:]
    p2 = segment[1,:]
    u = p2 - p1
    #print(u)
    norm_u = np.linalg.norm(u)
    
    p = np.zeros_like(q)
    dist = np.zeros(len(q))
    
    for i in range(len(q)):
        v = q[i] - p1
        check = np.dot(u, v) / (norm_u ** 2)
        
        if 0 <= check <= 1:
            p[i,:] = p1 + check * u
        elif check < 0:
            p[i,:] = p1
        else:
            p[i,:] = p2
        
        dist[i] = np.linalg.norm(q[i] - p[i])
    
    return p, dist




def compute_distance2fault2D(fault,stations,logic_sbnd,Merc,hemisphere='N'):
    """
    Calculate the closest point of the segment to point(s) in q

    Parameters:
    segment (numpy.ndarray(2,2)): Array containing the coordinates (easting, northing) of the end-points of the fault segment
    q (numpy.ndarray(N,2)): Array of coordinates (easting, northing) of points of interest (e.g., barycenters, stations).

    Returns:
    p (numpy.ndarray (N,2)): Array containing the coordinates (easting, northing) of the closest point of the segment to each point in q
    dist (numpy.ndarray (N)): Array of distance between the closest point of the segment to each point in q
    """
    #Adjust longitudes
    #fault[:,0][fault[:,0]>180]-=360  

    fault_UTM=np.zeros_like(fault)
    #fault_UTM[:,0],fault_UTM[:,1]=mesh.ll2xy(fault[:,1], fault[:,0], zone=Merc, hemisphere=hemisphere)
    fault_UTM[:,0],fault_UTM[:,1]=mesh.from_latlon(fault[:,1], fault[:,0], zone=Merc)
    fault_UTM[:,2]=fault[:,2]
    stations_UTM=np.zeros_like(stations)
    #stations_UTM[:,0],stations_UTM[:,1]=mesh.ll2xy(stations[:,1], stations[:,0], zone=Merc, hemisphere=hemisphere)
    stations_UTM[:,0],stations_UTM[:,1]=mesh.from_latlon(stations[:,1], stations[:,0], zone=Merc)
    stations_UTM[:,2]=stations[:,2]
    
    p0_aux=np.zeros([len(stations_UTM),len(fault_UTM)-1,2])
    dist_aux=np.zeros([len(stations_UTM),len(fault_UTM)-1])
    for j in range(len(fault_UTM)-1):
        p0_aux[:,j], dist_aux[:,j] = find_closest_point_segment(fault_UTM[[j,j+1],:][:,[0,1]], stations_UTM[:,[0,1]])
    
    distanceJB=np.min(dist_aux,axis=1)*1.e-3
    index_p0=np.argmin(dist_aux,axis=1)
    depth_fault=0.5*(fault_UTM[index_p0,2]+fault_UTM[index_p0+1,2])*1.e-3
    distanceJB=np.sqrt(distanceJB**2+(depth_fault-stations_UTM[:,2]*1.e-3)**2) #distance of each barycenter to closest point at the fault boundary
    p0=p0_aux[np.arange(p0_aux.shape[0]),index_p0,:] #UTM coordinates of closest points
    
    #Check if station is indeed inside the boundary when boundary is specified
    if logic_sbnd:
        boundary_polygon = Polygon(fault)
        # Check if each point is inside the polygon using shapely's vectorized contains function
        inside_mask = contains(boundary_polygon, stations[:, 0], stations[:, 1])
        distanceJB[~inside_mask]=0

    return p0, distanceJB



def Rupture_area_cells(slab, i, j):
    
    EToE = slab.EToE.astype(int)
    Area_cells = slab.Area_cells
    AreaSL = slab.AreaSL
    WidthSL = slab.WidthSL
    index_magnitude = slab.index_magnitude
    barycenter = slab.barycenter
    Nodes = slab.nodes
    Cells = slab.cells
    Matrix_distance = slab.Matrix_distance
    shape = slab.shape
    Fact_Area = slab.Fact_Area

    fact_mu_z = np.ones(len(Area_cells))
    exponent_scaling = 1

    if slab.Stress_drop_logic:
        gamma1 = slab.gamma1[0][j]  #FIXED by ANTONIO SCALA 2024/07/17
        gamma2 = slab.gamma2[0][j]  #FIXED by ANTONIO SCALA 2024/07/17
        sum_gamma = gamma1 + gamma2
        exponent_scaling = 0.5 * (gamma2 / sum_gamma)
        fact_mu_z = slab.fact_mu_z[j][:] #corrected by ANTONIO SCALA

    Event = [{'cell': [], 'true': True, 'nodes4events': []} for _ in range(len(barycenter[i][j]))]

    for l, index in enumerate(barycenter[i][j]):
        
        Event[l]['true'] = True
        Area_event = 0
        ng = 0
        Event[l]['cell'].append(int(index))
        fact_LWA = fact_mu_z[index]  #corrected by ANTONIO SCALA
        Event[l]['nodes4events'] = Cells[Event[l]['cell'][ng]].tolist() #nodes of the cell with index_cell
        ng += 1
        k = 1
        #k_nodes = 4
        Area_event = Area_cells[index]
        get_out = False
        
        while not get_out: 
            
            for jj in range(3): #loop over adjacent cells of index_cell
                
                
                if EToE[Event[l]['cell'][ng - 1], jj] != 0:
                    
                    index_cell = np.where(Event[l]['cell'][0:k ] == EToE[Event[l]['cell'][ng - 1], jj]-1)[0]
                    
                    if len(index_cell) == 0:
                        
                        Event[l]['cell'].append(EToE[Event[l]['cell'][ng - 1], jj]-1)
                        
                        if slab.Sub_boundary_logic:
                            point = Point(slab.barycenters_all[Event[l]['cell'][k], 0], slab.barycenters_all[Event[l]['cell'][k], 1])
                            polygon = Polygon(slab.bnd_mesh)
                            in_polygon = polygon.contains(point)
                            #in_polygon = contains(polygon, slab.barycenters_all[Event[l]['cell'][k], 0], slab.barycenters_all[Event[l]['cell'][k], 1])
                        else:
                            in_polygon = True
                        
                        if in_polygon:
                            
                            Event[l]['nodes4events'].extend(Cells[Event[l]['cell'][k]].tolist())
                            
                            depth4event = Nodes[np.array(Event[l]['nodes4events'])-1, 2]
                            #max_depth = np.max(depth4event)
                            #min_depth = np.min(depth4event)
                            distance = 1.e-3 * Matrix_distance[Event[l]['nodes4events'][np.argmax(depth4event)]-1, Event[l]['nodes4events'][np.argmin(depth4event)]-1]
                            
                            if shape == 'Circle':
                                Area_event += Area_cells[Event[l]['cell'][k]]
                                #Area_event += Area_cells[EToE[Event[l]['cell'][ng - 1], jj]]
                                k += 1
                                #k_nodes += 3
                                
                                
                            elif shape == 'Rectangle':
                                
                                if distance < max([1.1 * np.sqrt(Fact_Area) * ((fact_LWA)**exponent_scaling) * WidthSL[index_magnitude[i], j], 50]):
                                    Area_event += Area_cells[Event[l]['cell'][k]]
                                    #Area_event += Area_cells[EToE[Event[l]['cell'][ng - 1], jj]]
                                    k += 1
                                    #k_nodes += 3
                                else:
                                    Event[l]['cell'].pop()
                                    Event[l]['nodes4events']=Event[l]['nodes4events'][:-3]
                        else:
                            Event[l]['cell'].pop()

            
                if Area_event > Fact_Area * np.sqrt(fact_LWA) * AreaSL[index_magnitude[i], j]:
                    get_out = True
                    
                elif jj == 2:
                    ng += 1
                    
                    if ng - 1 > k or ng - 1 > len(slab.barycenters_all)-1 or ng - 1 > len(Event[l]['cell'])-1:
                        Event[l]['true'] = False
                        get_out = True
                        
    #print(Area_event)                  
    return Event



