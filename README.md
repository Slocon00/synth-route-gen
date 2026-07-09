# Synthetic route generation

The code in this repository is an extension of an existing project developed by the KDD Lab at ISTI-CNR (Pisa). The goal of the project is to generate and clone Individual Mobility Networks from a source area, where trajectory data is available, to a destination area, where such data is potentially missing. This extension focuses on route generation to connect the IMN locations mapped onto the destination area, implementing a more sophisticated approach to learn user preferences and compute user-customized routes (instead of relying on a simple, unrealistic shortest path computation).

Part of the code uses the Valhalla routing engine for map-matching, which requires a running instance to generate tile extracts from OpenStreetMap; refer to the [Valhalla documentation](https://valhalla.readthedocs.io/en/latest/) for more information on how to set up and run Valhalla.

## Requirements
```
Python==3.10.9
Flask==3.1.0
folium==0.18.0
geohash2==1.1
geopandas==1.0.1
geopy==2.4.1
networkx==3.4.2
numpy==2.1.3
osmnx==2.0.0
pandas==2.2.3
pytz==2024.2
rasterio==1.4.2
Requests==2.32.3
scikit_learn==1.5.2
scipy==1.14.1
Shapely==2.0.6
tqdm==4.67.1
pyvalhalla==3.2.0
```

## Running the code
Below is a brief description of the main functions and how to run them. Refer to the `example.ipynb` notebook for interactive execution of the code.

### Loading the data
Trajectory data is assumed to be GPS traces, in tabular format, with columns `userid`, `timestamp`, `latitude`, and `longitude`.

### Building the IMN


### Data exploration