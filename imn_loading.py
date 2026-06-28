import gzip
import json
import numpy as np
import pandas as pd
from math import radians, cos, sin, asin, sqrt
import osmnx as ox
from itertools import chain

def read_imn_old_format(filepath, period="03-04"):
    user_imns = {}
    imn_filedata = gzip.GzipFile(filepath, 'r')
    for row in imn_filedata:
        if len(row) <= 1:
            print('(skipping an empty IMN...)')
            continue
        
        user_obj = json.loads(row)
        uid = user_obj['uid']
        if period in user_obj:
            imn = user_obj[period]
            if imn is not None:
                # Extract key information for mobility cloning
                # At the present, only locations and real trips info
                imn_cloning = { 
                    'locations': { loc: {
                        'coordinates': imn['location_prototype'][loc],
                        'frequency': imn['location_features'][loc]['loc_support']
                        } for loc in imn['location_prototype'].keys() },
                    'trips': [(str(imn['traj_location_from_to'][ii][0]),  # origin location
                               str(imn['traj_location_from_to'][ii][1]),  # destination
                               imn['tid_se_times'][ii][0],           # start time
                               imn['tid_se_times'][ii][1],           # end time
                              ) for ii,_ in imn['traj_location_from_to'].items()],
                    'home': '0',  # simple frequency-based. TODO: use timeslot_count
                    'work': '1'
                }                
                user_imns[uid] = imn_cloning
    return user_imns


def get_trajectories_from_imns(imns: dict) -> dict:
    """
    Recovers the segmented trajectories from Individual Mobility Networks
    (IMN) in a separate dictionary.

    Parameters
    ----------
    imns : dict
        The IMN dictionary, as returned by read_imn(), containing multiple user
        IMNs.

    Returns
    -------
    segmented_trajs : dict
        A dictionary mapping each uid to a pandas DataFrame with columns
        'id', 'longitude', 'latitude', 'timestamp', representing the user's
        trajectories.
    """
    segmented_trajs = {
        np.int64(uid): pd.DataFrame(
                        data=np.array([
                            (traj['id'], p[0], p[1], p[2])
                            for traj in imn_data['trajectories'].values()
                            for p in traj['object']
                        ]),
                        columns=['id', 'longitude', 'latitude', 'timestamp']
                    ).astype({'id': int, 'timestamp': int})
        for uid, imn_data in imns.items()
    }
    return segmented_trajs


def read_imn(filepath):
    user_imns = {}
    imn_filedata = gzip.GzipFile(filepath, 'r')
    for row in imn_filedata:
        if len(row) <= 1:
            print('(skipping an empty IMN...)')
            continue
        
        imn = json.loads(row)
        uid = imn['uid']
        del imn['uid']
        if len(imn.keys()) != 0:
            # Extract key information for mobility cloning
            # At the present, only locations and real trips info
            imn_cloning = { 
                'locations': { int(loc): {
                    'coordinates': imn['location_prototype'][loc],
                    'frequency': imn['location_features'][loc]['loc_support']
                    } for loc in imn['location_prototype'].keys() },
                'trips': [(str(imn['traj_location_from_to'][ii][0]),  # origin location
                           str(imn['traj_location_from_to'][ii][1]),  # destination
                           imn['tid_se_times'][ii][0],           # start time
                           imn['tid_se_times'][ii][1],           # end time
                          ) for ii,_ in imn['traj_location_from_to'].items()],
                'home': '0',  # simple frequency-based. TODO: use timeslot_count
                'work': '1',
                'trajectories': imn['trajectories'],
                'o_d_trajectories': {
                    eval(k) : v for k, v in imn['location_traj_from_to'].items()
                },
                'trajectories_o_d': {
                    int(k): v for k, v in imn['traj_location_from_to'].items()
                },
            }                
            user_imns[uid] = imn_cloning
    return user_imns
