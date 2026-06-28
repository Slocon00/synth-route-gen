# %%
import numpy as np
import osmnx as ox
import networkx as nx
from shapely.geometry import box
from rtree import index
import gzip
import json


input_filename = "data/emilia_imns_UIDis00.json.gz"
output_filename = "data/emilia_poi_UIDis00.json.gz"
buffer_poi_meters = 200

# Open and read the gzipped JSON file
data = []
with gzip.open(input_filename, 'rt', encoding='utf-8') as f:
    lcount = 0
    for line in f:
        data.append(json.loads(line))
        lcount += 1
        print(f"reading IMN n.{lcount}")
        #break

# %%
loc_bbox = [999, 999, -999, -999]
for mydata in data:
    xlist = [xy[0] for _, xy in mydata['location_prototype'].items()]
    ylist = [xy[1] for _, xy in mydata['location_prototype'].items()]
    loc_bbox[0] = min(loc_bbox[0], min(xlist))
    loc_bbox[1] = min(loc_bbox[1], min(ylist))
    loc_bbox[2] = max(loc_bbox[2], max(xlist))
    loc_bbox[3] = max(loc_bbox[3], max(ylist))
loc_bbox = [loc_bbox[0]-0.001, loc_bbox[1]-0.001, loc_bbox[2]+0.001, loc_bbox[3]+0.001]  # data bbox + ~100m buffer
print(f"BBox of data: {loc_bbox}")


# %%
# Compute a meters_to_degree conversion factor to generate bbox guaranteed to contain all the data
# Looks at the 4 corners of the bbox and takes the worst case
xa, ya, xb, yb = loc_bbox
meters_per_cent_degree = []
for pt in [(xa,ya), (xa,yb), (xb,ya), (xb,yb)]:
    for delta in [(-0.001, 0), (0.001, 0), (0, -0.001), (0, 0.001)]:
        pt2 = (pt[0]+delta[0], pt[1]+delta[1])
        meters_per_cent_degree.append(ox.distance.great_circle(pt[1], pt[0], pt2[1], pt2[0]))
meters_per_cent_degree = np.min(meters_per_cent_degree)
meters_to_degree = 1/(meters_per_cent_degree/0.001) # factor to obtain meters from degrees


# %%
# Retrieve all amenities from OSMnx for the new bounding box
print("Retrieving amenities from OSMnx...")
tags = {'amenity': True}
gdf_amenities = ox.features.features_from_bbox(bbox=loc_bbox, tags=tags)
print(f" - Retrieved {len(gdf_amenities)} amenities:")

# %%
# Define classes for amenities
# Missing ones are considered "others"
amenity_classes = {
    'transportation': [
        'parking', 'taxi', 'fuel', 'parking_entrance', 'parking_exit', 'parking_access', 'car_rental', 
        'bicycle_rental', 'bus_station', 'bicycle_parking', 'motorcycle_parking', 'ferry_terminal', 
        'kick-scooter_parking', 'motorcycle_rental', 'taxi_rank'
    ],
    'healthcare': [
        'pharmacy', 'clinic', 'hospital', 'dentist', 'doctors', 'veterinary', 'nursing_home', 
        'health_post'
    ],
    'public_services': [
        'post_office', 'police', 'fire_station', 'townhall', 'courthouse', 'waste_disposal', 
        'waste_basket', 'social_facility', 'shelter', 'community_centre', 'public_bookcase', 
        'animal_shelter', 'events_venue', 'prison', 'archive', 'post_depot', 'mortuary', 
        'public_building', 'crematorium', 'payment_centre', 'meeting_room', 'reception_desk', 
        'reception_desk;lost_property_office', 'lost_property_office', 'group_home', 'dormitory'
    ],
    'finance': [
        'bank', 'atm', 'bureau_de_change', 'money_transfer'
    ],
    'food_and_drink': [
        'restaurant', 'cafe', 'fast_food', 'pub', 'bar', 'ice_cream', 'biergarten', 'food_court', 
        'canteen', 'restaurant;cafe'
    ],
    'utilities': [
        'telephone', 'recycling', 'drinking_water', 'toilets', 'vending_machine', 'charging_station', 
        'compressed_air', 'sanitary_dump_station', 'water_point', 'bicycle_repair_station', 
        'foot_shower', 'vacuum_cleaner', 'parcel_locker', 'fixme'
    ],
    'education': [
        'library', 'school', 'kindergarten', 'college', 'university', 'music_school', 'driving_school', 
        'childcare', 'research_institute', 'prep_school', 'language_school', 'training', 
        'dancing_school'
    ],
    'entertainment_and_recreation': [
        'cinema', 'theatre', 'nightclub', 'club', 'studio', 'events_venue', 'gambling', 'dojo', 
        'hookah_lounge', 'stripclub', 'concert_hall', 'planetarium', 'exhibition_centre', 
        'love_hotel', 'public_bath', 'watering_place', 'stage', 'auditorium'
    ],
    'shopping': [
        'marketplace', 'shop', 'mall', 'supermarket'
    ]
#    ,
#    'miscellaneous': [
#        'bench', 'fountain', 'place_of_worship', 'funeral_hall', 'car_wash', 'townhall', 'weighbridge', 'clock', 
#        'trade_union', 'coworking_space', 'printer', 'arts_centre', 'lavoir', 'shower', 'ticket_validator', 
#        'binoculars', 'loading_dock', 'grave_yard', 'luggage_locker', 'vehicle_inspection', 
#        'animal_boarding', 'left_luggage', 'waiting_room', 'dog_parking', 'office', 'sport_club', 
#        'monastery', 'boat_storage', 'press', 'factory', 'associaton', 'training', 'trolley_bay', 
#        'copyshop', 'dressing_room', 'table', 'smoking_area'
#    ]
}

# Create a reverse mapping from amenity type to class
amenity_to_class = {}
amenity_class_to_class_id = { class_name: id+1 for id, class_name in enumerate(amenity_classes.keys())  }
amenity_class_to_class_id['other'] = 0
class_id_to_amenity_class = { id: class_name for class_name, id in amenity_class_to_class_id.items() }
for class_name, amenities in amenity_classes.items():
    for amenity in amenities:
        amenity_to_class[amenity] = class_name

# %%
print("Loading amenities in R-tree index...")
# Create an empty R-tree index
idx = index.Index()

# Populate the R-tree index with the bounds of the amenities
index_n = 0
for _pos, amenity in gdf_amenities.iterrows():
    class_name = amenity_to_class.get(amenity['amenity'], 'other')
    class_id = amenity_class_to_class_id[class_name]
    idx.insert(index_n, amenity['geometry'].bounds, obj=class_id)
    index_n += 1

# %%
buffer_poi = buffer_poi_meters*meters_to_degree # meters

print("Computing POI frequencies...")
poi_classes = [class_id_to_amenity_class[i] for i in range(len(amenity_class_to_class_id))]
poi_dict_list = []
for mydata in data:
    # poi_dict associates each location to a vector of frequencies of amenities
    poi_dict = { 
        'uid': mydata['uid'], 
        'poi_classes': poi_classes,
        'poi_freq': { } 
        }
    for loc_id, xy in mydata['location_prototype'].items():
        x, y = xy
        x0, y0, x1, y1 = x-buffer_poi, y-buffer_poi, x+buffer_poi, y+buffer_poi 
        poi_vector = np.zeros(len(amenity_class_to_class_id))
        for item in list(idx.intersection((x0, y0, x1, y1), objects=True)):
            class_id = item.object
            # print(f"class_id: {class_id}, {item.object}, {item.bbox}")
            bbox = item.bbox
            if ox.distance.great_circle(y, x, bbox[1], bbox[0]) <= buffer_poi_meters:
                poi_vector[class_id] += 1
        poi_dict['poi_freq'][loc_id] = poi_vector.tolist()
    poi_dict_list.append(poi_dict)
    print(f" ({len(poi_dict_list)}/{len(data)}) uid: {mydata['uid']}, {len(poi_dict['poi_freq'])} locations processed")
    

# %%
print("Saving POI frequencies to file...")
with gzip.GzipFile(output_filename, 'w') as fout:
    for poi_dict in poi_dict_list:
        json_str = '%s\n' % json.dumps(poi_dict)
        json_bytes = json_str.encode('utf-8')
        fout.write(json_bytes)
