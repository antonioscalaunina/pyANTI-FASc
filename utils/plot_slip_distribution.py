import os
import glob
import json
from pathlib import Path
import time
from time import sleep
import folium
from branca.colormap import linear

# Define function to read JSON file
def read_config_json(file_path):
    with open(file_path, 'r') as file:
        return json.load(file)

def progress(percent=0, width=30):
    left = width * percent // 100
    right = width - left
    print('\r[', '#' * left, ' ' * right, ']',
          f' {percent: .0f}%', sep='', end='', flush=True)


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

def plot_slip_map(geojson_file, map_file, hypo):
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
    
    folium.Marker(
            location=[hypo[1], hypo[0]],
            popup="Epicenter",
            icon=folium.Icon(icon='star', color='black', prefix='fa')
        ).add_to(m)

    # Salva la mappa in un file HTML
    m.save(map_file)



# Clear screen
#os.system('cls' if os.name == 'nt' else 'clear')

# Read parameters from JSON file
with open('../config_files/Parameters/input.json') as fid:
    Param = json.load(fid)
Zone = Param['acronym']

# Define folder path
folder = '../output/'

# Check if application is 'PTF' and set variable accordingly
if Param['Configure']['application'] == 'PTF':
    hypo = Param['Event']['Hypo_LonLat']

# Define choices
values = ["event", "rigidity distribution", "magnitude", "scaling law"]

# Iterate through each choice
for ilev, value in enumerate(values, 1):
    print(f"Current folder is '{folder}'\n")
    choices = list(filter(lambda x: os.path.isdir(os.path.join(folder, x)), os.listdir(folder)))

    if len(choices) == 1:
        print(f"There is only one {value} directory\n")
        folder = os.path.join(folder+ choices[0]+ '/')
        print(folder)
        time.sleep(1.5)
    else:
        while True:
            print(f"Choose your {value} directory between:")
            for ndir, choice in enumerate(choices, 1):
                print(f"{ndir}. {choice}/")
            
            idir = input(f"Insert a number between 1 and {len(choices)}:\n\n")
            try:
                idir = int(idir)
                if 1 <= idir <= len(choices):
                    break
                else:
                    print(f"Please enter a number between 1 and {len(choices)}")
            except ValueError:
                print(f"Please enter a number between 1 and {len(choices)}")

        folder = os.path.join(folder+ choices[idir - 1]+ '/')

    # Clear screen
    # os.system('cls' if os.name == 'nt' else 'clear')

##############################################################

home=os.getcwd()
file_arr = glob.glob(folder+ 'Slip4H*.dat')
print ('Number of files: ', len(file_arr))
os.chdir(folder)
for i in range(len(file_arr)):
    progress(int(i*100/(len(file_arr)-1)))
    sleep(0.1)
    ascii_input_file = os.path.basename(file_arr[i])  # Assicurati di sostituire "input.txt" con il nome del tuo file ASCII
    geojson_output_file = (ascii_input_file[:-4]+'.json')  # Sostituisci con il nome del file GeoJSON che desideri creare
    map_file = (ascii_input_file[:-4]+'.html')	
    ascii_to_geojson(ascii_input_file, geojson_output_file)
    plot_slip_map(geojson_output_file, map_file, hypo)
os.chdir(home)
print("\n")
