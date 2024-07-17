"""
June 2024

"""
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
"""
import os
import glob
import json
from pathlib import Path
import time
from time import sleep
import folium
from branca.colormap import linear

"""



def get_min_dist_deg(slab):

    slab._check_attribute('EToE')
    dist=[]
    for i, element in enumerate(slab.EToE):
        lon,lat,_=slab.get_barycenter_coordinates([i])
        element=np.array([int(el-1) for el in element if el>0])
        lon_n, lat_n,_=slab.get_barycenter_coordinates(element)
        
    
        
        for i in range(len(lon_n)):
            dist.append(np.sqrt((lon_n[i]-lon)**2+(lat_n[i]-lat)**2))
    
    return max(dist)


def set_extent(lon,lat,ax,pad=1):

    lon_min=round(min(lon),1)-pad
    lon_max=round(max(lon),1)+pad
    lat_min=round(min(lat),1)-pad
    lat_max=round(max(lat),1)+pad
    ax.set_extent([lon_min, lon_max, lat_min, lat_max])

    return ax


def get_marker_size(slab,ax,fig):

    x1,x2=ax.get_xlim()
    y1,y2=ax.get_ylim()
    fig_width_in,fig_height_in=fig.get_size_inches()
    _, _, ax_width, ax_height = ax.get_position().bounds
    dpi = fig.dpi

    #size depending on axis limits
    ax_width_in=fig_width_in*ax_width
    x_points=ax_width_in*72
    ax_height_in=fig_height_in*ax_height
    y_points=ax_height_in*72
    x_deg=abs(x1-x2)
    y_deg=abs(y1-y2)
    
    dist_deg=get_min_dist_deg(slab)
    _15km_points=max(dist_deg*x_points/x_deg,dist_deg*y_points/y_deg)
    markersize_mesh=_15km_points**2

    markersize_single=max(ax_width_in, ax_height_in)/7
    

    return markersize_mesh, markersize_single

def get_scaled_fontsize(ax,fig):

    fig_width_in,fig_height_in=fig.get_size_inches()
    _, _, ax_width, ax_height = ax.get_position().bounds
    dpi = fig.dpi
    
    ax_width_in=fig_width_in*ax_width
    x_points=ax_width_in*72
    
    ax_height_in=fig_height_in*ax_height
    y_points=ax_height_in*72

    fontsize_axs=max(y_points,x_points)/40
    fontisize_cbar_ticks=max(y_points,x_points)/(0.7*50)
    fontisize_cbar_label=max(y_points,x_points)/(0.7*40)

    return fontsize_axs, fontisize_cbar_label, fontisize_cbar_ticks


def set_cbar(sc,ax,fig,label="",scale=None):

    _, _, ax_width, ax_height = ax.get_position().bounds
 
    cbar = plt.colorbar(sc, ax=ax, orientation='vertical', pad=0.05, aspect=30, shrink=0.7*ax_height)
    _,labelsize,labelsizeticks=get_scaled_fontsize(ax,fig)
    cbar.ax.tick_params(labelsize=labelsizeticks)
    cbar.set_label(label,size=labelsize)
    

    if scale is not None:
        
        ticks = cbar.get_ticks()
        scaled_tick_labels = [f'{tick * 10**scale:.2f}' for tick in ticks]
        # Set the new tick labels
        cbar.ax.set_yticklabels(scaled_tick_labels)
        cbar.ax.text(1.5, 1.05, rf'$\times 10^{{-{scale}}}$', transform=cbar.ax.transAxes, ha='center', va='bottom',size=labelsizeticks)
        



    return cbar



def set_dfbasemap(fig,ax):
    """
    Set and plot a default basemap in the Geoaxs ax

    """
    # Add coastlines and features
    ax.coastlines(resolution='10m',linewidth=0.7)
    ax.add_feature(cfeature.LAND, linestyle=':',color='gainsboro')
    ax.add_feature(cfeature.OCEAN, linestyle=':',color='#0085BD')
    
    # Add gridlines
    gl = ax.gridlines(draw_labels=True, alpha=0)
    gl.top_labels = False
    gl.right_labels = False
    

    fontsize,_,_=get_scaled_fontsize(ax,fig)
    gl.xlabel_style = {'size': fontsize}
    gl.ylabel_style = {'size': fontsize}
    
    return fig, ax

def get_slip(input_file):
    
    slip_a=[]
    with open(input_file, 'r') as f:
        next(f)  
        for line in f:
            data = line.strip().split()
            lon1, lat1, depth1, lon2, lat2, depth2, lon3, lat3, depth3, rake, slip = map(float, data)
            slip_a.append(slip)
    return slip_a


"""
# Define ascii to geojson
def ascii_to_geojson(input_file, output_file):
    features = []

    with open(input_file, 'r') as f:
        next(f)  # Salta la prima riga (header)
        for line in f:
            data = line.strip().split()
            lon1, lat1, depth1, lon2, lat2, depth2, lon3, lat3, depth3, rake, slip = map(float, data)

            # Costruzione del poligono GeoJSON
            polygon = {
                "type": "Feature",
                "properties": {
                    "rake": rake,
                    "slip": slip
                },
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [[
                        [lon1, lat1, depth1],
                        [lon2, lat2, depth2],
                        [lon3, lat3, depth3],
                        [lon1, lat1, depth1]  # Aggiungere nuovamente il primo punto per chiudere il poligono
                    ]]
                }
            }
            features.append(polygon)

    # Creazione dell'oggetto GeoJSON
    geojson_data = {
        "type": "FeatureCollection",
        "features": features
    }

    # Scrittura su file
    with open(output_file, 'w') as f:
        json.dump(geojson_data, f, indent=2)

def plot_slip_map(geojson_file, map_file):
    # Carica il file GeoJSON
    with open(geojson_file, 'r') as f:
        data = json.load(f)
    
    # Estrai le coordinate dal GeoJSON per determinare la posizione e lo zoom della mappa
    latitudes = []
    longitudes = []
    for feature in data['features']:
        coords = feature['geometry']['coordinates'][0]
        for lon, lat, _ in coords:
            latitudes.append(lat)
            longitudes.append(lon)

    # Determina la posizione centrale della mappa e il livello di zoom
    avg_lat = sum(latitudes) / len(latitudes)
    avg_lon = sum(longitudes) / len(longitudes)
    max_lat = max(latitudes)
    min_lat = min(latitudes)
    max_lon = max(longitudes)
    min_lon = min(longitudes)

    # Calcola il centro e lo zoom della mappa
    center_lat = (max_lat + min_lat) / 2
    center_lon = (max_lon + min_lon) / 2
    zoom = 10  # Puoi regolare lo zoom a tua discrezione

    # Imposta il centro e lo zoom della mappa
    #m.location = [center_lat, center_lon]
    #m.zoom_start = zoom

    # Inizializza la mappa Folium
    m = folium.Map(location=[center_lat, center_lon], zoom_start=7)

    # Definisci una scala di colori per lo slip
    cmap = linear.YlOrRd_09.scale(min(data['features'], key=lambda x: x['properties']['slip'])['properties']['slip'],
                                   max(data['features'], key=lambda x: x['properties']['slip'])['properties']['slip'])

    # Aggiungi i poligoni alla mappa
    for feature in data['features']:
        coords = feature['geometry']['coordinates'][0]  # Ottieni le coordinate del poligono
        slip = feature['properties']['slip']  # Ottieni lo slip
        depth = feature['geometry']['coordinates'][0][0][2]

        # Inversione della convenzione [longitudine, latitudine] in [latitudine, longitudine]
        coords = [[lat, lon] for lon, lat, _ in coords]

        # Crea il poligono e aggiungilo alla mappa con il colore corrispondente allo slip per il riempimento
        folium.Polygon(
            locations=coords,
            color='',  # Colore del bordo
            fill_color=cmap(slip),  # Colore di riempimento
            fill_opacity=0.7,  # Opacità del colore di riempimento
            tooltip=f"Slip: {slip}"
        ).add_to(m)

    # Aggiungi la legenda della scala di colori alla mappa
    cmap.caption = 'Slip (m)'
    cmap.add_to(m)

    # Salva la mappa in un file HTML
    m.save(map_file)
"""