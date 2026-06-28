import rasterio
from rasterio.mask import mask
import geopandas as gpd
import pandas as pd
import numpy as np
from shapely.geometry import Polygon
from shapely.geometry import Point
import osmnx as ox

def fetch_population_data_from_tiff(tiff_file, bbox):
    """
    Fetches population data from a local GeoTIFF file within a given bounding box.
    
    Parameters:
    tiff_file (str): Path to the GeoTIFF file containing population data
    bbox (tuple): Bounding box in the format (min_lon, min_lat, max_lon, max_lat)
    
    Returns:
    geopandas.GeoDataFrame: GeoDataFrame containing population data in grid format
    """
    global my_out_trans
    min_lon, min_lat, max_lon, max_lat = bbox

    # Open the GeoTIFF file
    with rasterio.open(tiff_file) as src:
        # Define the bounding box as a GeoJSON-like geometry
        bbox_geom = {
            "type": "Polygon",
            "coordinates": [[
                (min_lon, min_lat), (max_lon, min_lat), (max_lon, max_lat),
                (min_lon, max_lat), (min_lon, min_lat) ]] }
        
        # Use rasterio's mask function to crop the data to the bounding box
        out_image, out_transform = mask(src, [bbox_geom], crop=True)
        out_image = out_image[0]  # Assuming single band (population density data)

        # Extract cells with population data and create polygons
        population_data = []
        rows, cols = np.where(out_image > 0)  # Only consider cells with non-zero population
        for row, col in zip(rows, cols):
            population_density = out_image[row, col]
            x, y = out_transform * (col, row)
            x2,y2 = out_transform * (col+1, row+1)
            # Create a polygon for each cell
            cell_polygon = Polygon([
                (x, y), (x2 , y), (x2, y2),
                (x, y2), (x, y) ])
            population_data.append({"population_density": population_density, "geometry": cell_polygon})

    # Convert the data to a GeoDataFrame
    population_gdf = gpd.GeoDataFrame(population_data, geometry="geometry")
    population_gdf.set_crs(src.crs, inplace=True)

    return population_gdf
   


def add_intersecting_nodes(gdf, graph):
    """
    Extends a GeoDataFrame with a column containing lists of node IDs 
    from an OSMnx graph that intersect with each geometry in the GeoDataFrame.
    
    Parameters:
    gdf (GeoDataFrame): The GeoDataFrame with geometry to check intersections
    graph (networkx.MultiDiGraph): The OSMnx graph containing nodes
    
    Returns:
    GeoDataFrame: Extended GeoDataFrame with a new column 'intersecting_nodes'
    """
    # Convert graph nodes to a GeoDataFrame
    nodes = ox.graph_to_gdfs(graph, nodes=True, edges=False)
    
    # Ensure the CRS matches between the graph nodes and the input GeoDataFrame
    if gdf.crs != nodes.crs:
        nodes = nodes.to_crs(gdf.crs)

    # For each geometry in gdf, find intersecting nodes
    intersecting_node_ids = []
    for geom in gdf.geometry:
        # Check which nodes intersect the geometry
        intersecting_nodes = nodes[nodes.intersects(geom)]
        intersecting_node_ids.append(intersecting_nodes.index.tolist())  # Add node IDs as list

    # Add the list of intersecting node IDs as a new column in the GeoDataFrame
    gdf = gdf.copy()  # Avoid modifying the original GeoDataFrame
    gdf['intersecting_nodes'] = intersecting_node_ids
    
    return gdf




#def find_in_cumulative(df, p_r):
#    # Use numpy's searchsorted to find the position where p_r should be inserted
#    idx = np.searchsorted(df['cumulative_p'], p_r, side='left')
#    
#    # Check if idx is valid and within the dataframe range
#    if idx < len(df):
#        return df.iloc[idx]['intersecting_nodes']
#    else:
#        return None  # In case p_r is larger than all cumulative_p values
#
#
#def select_random_nodes(gdf_cumulative_p, n=10):
#    node_list = []
#    for pr in np.random.random(n):
#        candidates = find_in_cumulative(gdf_cumulative_p, pr)
#        node_list.append(candidates[np.random.randint(len(candidates))])
#    return node_list
    


def generate_cumulative_map(population_tiff_file, G):
    bbox = ox.graph_to_gdfs(G, nodes=True, edges=False).total_bounds
    population_gdf = fetch_population_data_from_tiff(population_tiff_file, bbox)
    gdf_extended = add_intersecting_nodes(population_gdf, G)
    gdf_extended['p'] = gdf_extended['population_density'] / gdf_extended.population_density.sum()
    gdf_extended['cumulative_p'] = gdf_extended['p'].cumsum()
    gdf_cumulative_p = gdf_extended[['cumulative_p','intersecting_nodes']]
    return gdf_cumulative_p

    
    
    
    

