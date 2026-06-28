import pandas as pd
import numpy as np
import osmnx as ox

import imn_generation
import sys
import os

# Read data points and rename the essential columns: id, coordinates and timestamp
# Notice: the IMN generation tool needs timestamps expressed as Epoch time, in seconds
# Notice 2: file format of Octo data:
#     ID_ANONYMOUS,DAY,HH24,LATITUDE,LONGITUDE,SPEED,HEADING,QUALITY,ID_PANELSESSION,DELTAPOS
#     400,2019-05-18,14:33:01,43647819,11465733,114,358,3,1,2076

if len(sys.argv) < 3:
    sys.exit("Error: Please provide both input and output file names as command line arguments.")

input_file = sys.argv[1]
output_file = sys.argv[2]

if not os.path.isfile(input_file):
    sys.exit(f"Error: The input file '{input_file}' does not exist.")
output_file = sys.argv[2]

points_df = pd.read_csv(input_file).rename(columns={
    'ID_ANONYMOUS':'id', 
    'LONGITUDE':'longitude', 
    'LATITUDE':'latitude'
})

points_df = points_df[(points_df['DELTAPOS'] > 0) | (points_df['ID_PANELSESSION'] == 0)] # remove points with no movement

points_df['timestamp'] = points_df['DAY'].astype(str) + ' ' + points_df['HH24'].astype(str)
points_df['timestamp'] = pd.to_datetime(points_df['timestamp']).values.astype(np.int64) // 10 ** 9
points_df['longitude'] = points_df['longitude']/10**6
points_df['latitude'] = points_df['latitude']/10**6
points_df = points_df[['id', 'latitude', 'longitude', 'timestamp']]

# Use the dataframe with points and store the generated IMNs in a zipped json file
imn_generation.main_from_code(points_df, output_file)



