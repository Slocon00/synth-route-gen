from math import radians, cos, sin, asin, sqrt
import osmnx as ox
import networkx as nx
import random
import pandas as pd
import numpy as np

# Cache: (source_node, destination_node) -> path as list of nodes
osm_paths_cache = { }


def haversine_distance(x1, y1, x2, y2):
    R = 6371000 #meters
    lat_rad1 = radians(y1)
    lon_rad1 = radians(x1)
    lat_rad2 = radians(y2)
    lon_rad2 = radians(x2)
    return 2*R * asin(sqrt(sin((lat_rad2-lat_rad1)/2)**2 + cos(lat_rad1)*cos(lat_rad2)*(sin((lon_rad2-lon_rad1)/2)**2)))


###################################################################

def retrieve_osm(center, radius=20000):
    # input: center in the form (lat, lon), radius in meters
    G = ox.graph.graph_from_point(center, dist=radius, network_type="drive", simplify=False)
    G = ox.routing.add_edge_speeds(G)
    G = ox.routing.add_edge_travel_times(G)

    return G
    

###################################################################

def shortest_path_straight_line(target_osm, from_osm, to_osm, start_time):
    from_x = target_osm.nodes[from_osm]['x']
    from_y = target_osm.nodes[from_osm]['y']
    to_x = target_osm.nodes[to_osm]['x']
    to_y = target_osm.nodes[to_osm]['y']
    dist = haversine_distance(from_x, from_y, to_x, to_y)
    travel_time = int( dist / (50 / 3.6) )  # 50 km/h as m/s, travel time in secs
    n_samples = travel_time//10 + 1
    traj = []
    for i in range(n_samples):
        traj.append([start_time + i*10,
                     from_x + (to_x-from_x)*i/n_samples, 
                     from_y + (to_y-from_y)*i/n_samples
                    ])
    traj.append([start_time + travel_time, to_x, to_y])
    return traj, travel_time


###################################################################

def create_trajectory(G, node_origin, node_destination, start_time, weight="travel_time", slow_factor=1, verbose=False, use_cache=True):
    """
    Creates a trajectory between two nodes in a road network graph.
    
    Parameters:
        G (networkx.Graph): A road network graph.
        node_origin (int): The origin node ID.
        node_destination (int): The destination node ID.
        weight (str): The edge attribute to consider as the weight for computing the shortest path (default is "travel_time").
        slow_factor (float): A factor to adjust the travel time estimation (default is 1).
        
    Returns:
        tuple: A tuple containing:
            - trajectory (list): A list of tuples, each containing GPS coordinates and associated time for a point along the trajectory.
            - total_time (float): The total travel time of the trajectory.
            - road_segments_osmid (list): A list of road segment IDs along the trajectory.
    """
    global osm_paths_cache
    
    if verbose:
        print(f"Computing path: {node_origin} --> {node_destination}...")
        
    # If self-loop, assume a non-movement of 1 minute
    if node_origin == node_destination:
        return [(G.nodes[node_origin]["y"], G.nodes[node_origin]["x"], start_time), 
                (G.nodes[node_destination]["y"], G.nodes[node_destination]["x"], start_time + 60)], 60, []
    
    # Compute the shortest path between the two nodes according to a given weight
    if use_cache:
        if (node_origin,node_destination) in osm_paths_cache:
            route = osm_paths_cache[(node_origin, node_destination)]
        else:
            try:
                route = nx.shortest_path(G, node_origin, node_destination, weight=weight)
            except nx.NetworkXNoPath:
                route = None
            osm_paths_cache[(node_origin, node_destination)] = route

    else:
#        route = nx.shortest_path(G, node_origin, node_destination, weight=weight)
        try:
            route = nx.shortest_path(G, node_origin, node_destination, weight=weight)
        except nx.NetworkXNoPath:
            route = None
   
    # Routing failure !
    if route == None:
        return None, None, None
    
    
    # Retrieve information about the route
    gdf_info_routes = ox.routing.route_to_gdf(G, route)
        
    # Get the road segment IDs
    road_segments_osmid = list(gdf_info_routes[["osmid"]].index)
    
    # Extract sequence of GPS points
    gps_points = [(G.nodes[node]["y"], G.nodes[node]["x"]) for node in route]
    
    # Guess the travel time for each edge and add the temporal dimension to the GPS points
    travel_time_edges = gdf_info_routes["travel_time"].values
    travel_time_guess_edges = travel_time_edges * slow_factor
    list_seconds = [0] + list(travel_time_guess_edges.cumsum())
    
    # Combine GPS points with associated time
    trajectory = [(gps[0], gps[1], time + start_time) for gps, time in zip(gps_points, list_seconds)]
    
    return trajectory, list_seconds[-1], road_segments_osmid


###################################################################

def find_in_cumulative(df, p_r):
    # Used to sample random nodes following population density
    
    # Use numpy's searchsorted to find the position where p_r should be inserted
    idx = np.searchsorted(df['cumulative_p'], p_r, side='left')
    
    # Check if idx is valid and within the dataframe range
    if idx < len(df):
        return df.iloc[idx]['intersecting_nodes']
    else:
        return None  # In case p_r is larger than all cumulative_p values


###################################################################

def select_random_nodes(gdf_cumulative_p, n=10):
    # Used to sample random nodes following population density
    node_list = []
    for pr in np.random.random(n):
        candidates = find_in_cumulative(gdf_cumulative_p, pr)
        node_list.append(candidates[np.random.randint(len(candidates))])
    return node_list



###################################################################


def map_imn_to_osm(imn, target_osm, home_osm=None, work_osm=None, n_trials=10, gdf_cumulative_p=None):
    # mapping locations from imn to target_osm
    
    # Special case: home and work
    # home_osm and work_osm, if None, are obtained from "imn"
    best_dist = 999999999
    nodes = list(target_osm.nodes())
    if gdf_cumulative_p == None:
        gdf_cumulative_p = pd.DataFrame({ "cumulative_p": [1.0], "intersecting_nodes": [ nodes ]})
    imn_dist = haversine_distance(imn['locations'][imn['home']]['coordinates'][0],
                                 imn['locations'][imn['home']]['coordinates'][1],
                                 imn['locations'][imn['work']]['coordinates'][0],
                                 imn['locations'][imn['work']]['coordinates'][1])
    for i in range(n_trials):
        tmp_home_osm = home_osm
        if home_osm == None:
            tmp_home_osm = select_random_nodes(gdf_cumulative_p, n=1)[0]  # nodes[random.randint(0,len(nodes)-1)]
        tmp_work_osm = work_osm
        if work_osm == None:
            tmp_work_osm = select_random_nodes(gdf_cumulative_p, n=1)[0]  # nodes[random.randint(0,len(nodes)-1)]
        dist = haversine_distance(target_osm.nodes[tmp_home_osm]['x'],
                                 target_osm.nodes[tmp_home_osm]['y'],
                                 target_osm.nodes[tmp_work_osm]['x'],
                                 target_osm.nodes[tmp_work_osm]['y'])

        if abs(dist-imn_dist) < abs(best_dist-imn_dist):
            best_home = tmp_home_osm
            best_work = tmp_work_osm
            best_dist = dist

    map_loc = { imn['home']: best_home, imn['work']: best_work }
    SE = (best_dist/1000 - imn_dist/1000)**2  # Sum of squared errors
    
    # General case
    for loc in imn['locations'].keys():  # Notice: they are sorted by loc_id, thus frequency
        if loc in map_loc:
            continue
        best_dist = 99999999999
        for i in range(n_trials):
            tmp_osm = select_random_nodes(gdf_cumulative_p, n=1)[0]  # nodes[random.randint(0,len(nodes)-1)]
            dist = 0  # contains the sum of squared distance errors
            for imn_l, osm_l in map_loc.items():  # compare against the nodes already mapped
                dist_osm = haversine_distance(target_osm.nodes[tmp_osm]['x'],
                                              target_osm.nodes[tmp_osm]['y'],
                                              target_osm.nodes[osm_l]['x'],
                                              target_osm.nodes[osm_l]['y'])
                dist_imn = haversine_distance(imn['locations'][loc]['coordinates'][0],
                                              imn['locations'][loc]['coordinates'][1],
                                              imn['locations'][imn_l]['coordinates'][0],
                                              imn['locations'][imn_l]['coordinates'][1])
                dist += (dist_osm/1000 - dist_imn/1000)**2

            if dist < best_dist:
                best_osm = tmp_osm
                best_dist = dist
                best_i = i
        map_loc[loc] = best_osm
        SE += best_dist

    if len(map_loc) > 1:
        rmse = sqrt(SE/(len(map_loc)*(len(map_loc)-1)/2))
    else:
        rmse = 0

    return map_loc, rmse


###################################################################
def init_osm_paths_cache(imn, G, map_loc, use_prefetch=False):
    global osm_paths_cache
    
    if use_prefetch:
        for i in range(len(osm_paths_cache)-500000):  # Cap cache size (plus the new ones to add)
            osm_paths_cache.popitem()
        map_loc_set = [map_loc[loc] for loc in imn['locations'].keys()]
        for loc in imn['locations'].keys():
            if imn['locations'][loc]['frequency'] > 10:
                node_origin = map_loc[loc]
                routes = nx.single_source_dijkstra_path(G, source=node_origin, cutoff=None, weight="travel_time")
                for dest, path in routes.items():
                    if dest in map_loc_set:
                        osm_paths_cache[(node_origin,dest)] = path
    else:
        osm_paths_cache = {}


###################################################################


def simulate_trips(imn, target_osm, 
                   start_sym=None, end_sym=None, 
                   home_osm=None, work_osm=None, 
                   n_trials=10, stay_time_error=0.2,
                   gdf_cumulative_p=None,
                   use_cache=True,
                   use_prefetch=False):
    # Maps the trips in "imn" between timestamp "start_sym" and "end_sym"

    # randomly map locations
    map_loc, rmse = map_imn_to_osm(imn, target_osm, 
                             home_osm=home_osm, work_osm=work_osm, 
                             n_trials=n_trials,
                             gdf_cumulative_p=gdf_cumulative_p)
    
    # scan imn trips, selecting those within [start_sym, end_sym]
    if start_sym == None:
        start_sym = min([r[2] for r in imn['trips']]) - 1
    if end_sym == None:
        end_sym = max([r[2] for r in imn['trips']]) + 1
    
    sym_time = None
    prev_end = None
    
    trajectory = []
    osm_segments_usage = {} 
    if use_cache:
        init_osm_paths_cache(imn, target_osm, map_loc, use_prefetch=use_prefetch)  # each user starts from a scratch cache
    for l_from, l_to, st, end in imn['trips']:
        if l_from == l_to:
            continue
        if (st < start_sym) or (st > end_sym):
            continue
        if sym_time == None:
            sym_time = st
            prev_end = st
        # Use IMN stay time with randomization
        sym_stay = (st - prev_end) + random.randint(-int((st - prev_end)*stay_time_error), 
                                                    int((st - prev_end)*stay_time_error))
        sym_time += sym_stay
        prev_end = end
        
        # Shortest path with its duration
        shortest_p, duration_p, osm_segments = create_trajectory(target_osm, map_loc[l_from], map_loc[l_to], sym_time, use_cache=use_cache)
        
        # If shortest path fails, all the IMN mapping fails
        if shortest_p == None:
            return None, None, None, None
        
        trajectory.extend(shortest_p)  # TODO: resample trajectory to 1 Hz by linear interpolation
        sym_time += duration_p
        for oss in osm_segments:
            osm_segments_usage[oss] = osm_segments_usage.get(oss,0) + 1

    return trajectory, osm_segments_usage, map_loc, rmse
    
