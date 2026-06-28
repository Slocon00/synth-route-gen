import math
import pytz
import datetime

from collections import defaultdict

R_EARTH = 6371000


def dist2angle(dist):

    return dist * 180.0 / math.pi / R_EARTH


def get_ordered_history(imh):
    history_order_dict = dict()
    for tid in imh['trajectories']:
        ts = datetime.datetime.fromtimestamp(imh['trajectories'][tid].start_point()[2] / 1000.0)
        history_order_dict[tid] = ts
    return history_order_dict


def default_to_regular(d):
    if isinstance(d, defaultdict):
        d = {k: default_to_regular(v) for k, v in d.items()}
    return d


rome_params = {
    'input_table': 'agenda.rome1000',
    'min_lat': 41.24,
    'min_lon': 11.59,
    'tzoffset': int(datetime.datetime.now(pytz.timezone('Europe/Rome')).strftime('%z')[:3]),
    'traintest_date': datetime.datetime.strptime('2015-05-03', '%Y-%m-%d'),
}

london_params = {
    'input_table': 'agenda.london1000',
    'min_lat': 51.15,
    'min_lon': -0.89,
    'tzoffset': int(datetime.datetime.now(pytz.timezone('Europe/London')).strftime('%z')[:3]),
    'traintest_date': datetime.datetime.strptime('2015-05-03', '%Y-%m-%d'),
}

boston_params = {
    'input_table': 'agenda.boston1000',
    'min_lat': 40.91,
    'min_lon': -73.98,
    'tzoffset': int(datetime.datetime.now(pytz.timezone('US/Eastern')).strftime('%z')[:3]),
    'traintest_date': datetime.datetime.strptime('2015-05-03', '%Y-%m-%d'),
}

beijing_params = {
    'input_table': 'agenda.beijing',
    'min_lat': 38.945889,
    'min_lon': 114.698094,
    'tzoffset': int(datetime.datetime.now(pytz.timezone('Asia/Shanghai')).strftime('%z')[:3]),
    'traintest_date': datetime.datetime.strptime('2008-05-04', '%Y-%m-%d'),
}
