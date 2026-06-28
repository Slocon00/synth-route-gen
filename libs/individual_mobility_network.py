from tosca import *
from util import *
from mobility_distance_functions import trajectory_distance

import networkx as nx
import pandas as pd
import warnings
warnings.simplefilter("ignore")
import time


__author__ = 'Riccardo Guidotti'


def get_points_trajfromto(imh):

    #Added By Omid
    trajectories = imh['trajectories']
    #trajectories = imh

    points = dict()
    traj_from_to = dict()

    for tid, traj in trajectories.items():

        lon_from = float(traj.start_point()[0])
        lat_from = float(traj.start_point()[1])
        time_from = int(traj.start_point()[2])

        lon_to = float(traj.end_point()[0])
        lat_to = float(traj.end_point()[1])
        time_to = int(traj.end_point()[2])

        pid_start_point = len(points)
        points[pid_start_point] = [lon_from, lat_from, time_from, 'f', tid]

        pid_end_point = len(points)
        points[pid_end_point] = [lon_to, lat_to, time_to, 't', tid]

        traj_from_to[tid] = [pid_start_point, pid_end_point]

    return points, traj_from_to


def radius_of_gyration(points, centre_of_mass, dist):
    rog = 0
    for p in points:
        rog += dist(p, centre_of_mass)
    rog = 1.0*rog/len(points)
    return rog


def entropy(x, classes=None):
    if len(x) == 1:
        return 0.0
    val_entropy = 0
    n = np.sum(x)
    for freq in x:
        if freq == 0:
            continue
        p = 1.0 * freq / n
        val_entropy -= p * np.log2(p)
    if classes is not None and classes:
        val_entropy /= np.log2(classes)
    return val_entropy


def normalize_dict(x):

    max_val = np.max(list(x.values()))

    norm_dict = dict()
    for k in x:
        norm_dict[k] = 1.0 * x[k] / max_val

    return norm_dict

# def get_convex_hull(points):
#     p_list = list()
#
#     for p in points:
#         p_list.append((p[0], p[1]))
#
#     if len(p_list) < 3:
#         geos = GeoSeries([Point(p_list[0])])
#     else:
#         geos = GeoSeries([Polygon(p_list)])
#
#     gdf = gpd.GeoDataFrame({'geometry': geos.convex_hull})
#
#     convex_hull_shape = json.loads(gdf.to_json())['features'][0]['geometry']
#
#     return convex_hull_shape


def locations_detection(points, min_dist=50.0, nrun=5):

    # npoints = len(points)

    spatial_points = list()
    for p in list(points.values()):
        if p[0] in [np.nan, np.inf] or p[1] in [np.nan, np.inf]:
            continue
        spatial_points.append(p[0:2])

    if len(spatial_points) == 0:
        return None

    centers_min, centers_max = get_min_max(spatial_points)

    cluster_res = dict()
    cuts = dict()
    for runid in range(0, nrun):
        try:
            tosca = Tosca(kmin=centers_min, kmax=centers_max, xmeans_df=spherical_distances,
                          singlelinkage_df=spherical_distance, is_outlier=thompson_test,
                          min_dist=min_dist, verbose=False)
            tosca.fit(np.asarray(spatial_points))
            cluster_res[tosca.k_] = tosca.cluster_centers_
            cuts[tosca.k_] = tosca.cut_dist_
        except ValueError:
            pass

    if len(cluster_res) == 0:
        return None

    index = np.min(list(cluster_res.keys()))
    centers = cluster_res[index]
    loc_tosca_cut = cuts[index]

    # calculate distances between points and medoids
    distances = spherical_distances(np.asarray(spatial_points), np.asarray(centers))

    # calculates labels according to minimum distance
    labels = np.argmin(distances, axis=1)

    # build clusters according to labels and assign point to point identifier
    location_points = defaultdict(list)
    location_prototype = dict()
    for pid, lid in enumerate(labels):
        location_points[lid].append(pid)
        location_prototype[lid] = list(centers[lid])

    # rename locations from bigger to smaller
    pid_lid = dict()
    location_support = dict()
    location_sorted = sorted(location_points.keys(), key=lambda x: len(location_points[x]), reverse=True)
    new_location_points = defaultdict(list)
    new_location_prototype = defaultdict(list)
    for new_lid, lid in enumerate(location_sorted):
        new_location_points[new_lid] = location_points[lid]
        new_location_prototype[new_lid] = location_prototype[lid]
        location_support[new_lid] = len(location_points[lid])
        for pid in location_points[lid]:
            pid_lid[pid] = new_lid

    location_points = new_location_points
    location_prototype = new_location_prototype

    # statistical information for users analysis
    cm = np.mean(spatial_points, axis=0)
    rg = radius_of_gyration(spatial_points, cm, spherical_distance)
    en = entropy(list(location_support.values()), classes=len(location_support))

    res = {
        'location_points': location_points,
        'location_prototype': location_prototype,
        'pid_lid': pid_lid,
        'rg': rg,
        'entropy': en,
        'loc_tosca_cut': loc_tosca_cut,
    }

    return res


def movements_detection(pid_lid, traj_from_to, imh):

    traj_from_to_loc = dict()
    loc_from_to_traj = defaultdict(list)
    loc_nextlocs = defaultdict(lambda: defaultdict(int))

    for tid, from_to in traj_from_to.items():
        loc_from = pid_lid[from_to[0]]
        loc_to = pid_lid[from_to[1]]
        traj_from_to_loc[tid] = [loc_from, loc_to]
        loc_from_to_traj[(loc_from, loc_to)].append(tid)
        loc_nextlocs[loc_from][loc_to] += 1

    G = nx.DiGraph()
    for loc_from in loc_nextlocs:
        for loc_to in loc_nextlocs[loc_from]:
            nbr_traj = loc_nextlocs[loc_from][loc_to]
            G.add_edge(loc_from, loc_to, weight=nbr_traj)

    movement_traj = dict()
    lft_mid = dict()
    for mid, lft in enumerate(loc_from_to_traj):
        movement_traj[mid] = [lft, loc_from_to_traj[lft]]
        lft_mid[lft] = mid

    #print(loc_from_to_traj)
    #print(movement_traj)
    trajectories = imh['trajectories']
    movement_prototype = dict()

    for mid in movement_traj:
        traj_in_movement = movement_traj[mid][1]

        if len(traj_in_movement) > 2:
            prototype = None
            min_dist = float('inf')
            for tid1 in traj_in_movement:
                tot_dist = 0.0
                traj1 = trajectories[tid1]
                for tid2 in traj_in_movement:
                    traj2 = trajectories[tid2]
                    dist = trajectory_distance(traj1, traj2)
                    tot_dist += dist
                if tot_dist < min_dist:
                    min_dist = tot_dist
                    prototype = traj1
            movement_prototype[mid] = prototype
        else:
            movement_prototype[mid] = trajectories[traj_in_movement[0]]

    res = {
        'movement_traj': movement_traj,
        'movement_prototype': movement_prototype,
        'loc_nextlocs': loc_nextlocs,
        'traj_from_to_loc': traj_from_to_loc,
        'lft_mid': lft_mid,
        'graph': G,
    }

    return res


def get_location_features(points_in_loc, traj_from_to_loc, location_prototype, imh):

    trajectories = imh['trajectories']
    sorted_points = sorted(points_in_loc, key=lambda x: points_in_loc[x][2])

    staytime_dist = defaultdict(int)
    nextloc_count = defaultdict(int)
    nextloc_dist = defaultdict(lambda: defaultdict(int))
    staytime_durations = list()

    for i in range(0, len(sorted_points)-1):
        pid1 = sorted_points[i]
        pid2 = sorted_points[i+1]

        arriving_leaving = points_in_loc[pid1][3]

        if arriving_leaving == 't':

            ts1 = datetime.datetime.fromtimestamp(points_in_loc[pid1][2])
            ts2 = datetime.datetime.fromtimestamp(points_in_loc[pid2][2])

            at = ts1.replace(second=0, microsecond=0)
            lt = ts2.replace(second=0, microsecond=0)

            midnight_at = at.replace(hour=0, minute=0)
            midnight_lt = lt.replace(hour=0, minute=0)

            at_sec = int((at - midnight_at).total_seconds())
            lt_sec = int((lt - midnight_lt).total_seconds())

            if at_sec <= lt_sec:
                for minute in range(at_sec, lt_sec + 60, 60):
                    dt_minute = datetime.time(hour=int(minute/3600), minute=int((minute % 3600)/60))
                    staytime_dist[dt_minute] += 1
            elif at_sec > lt_sec:
                for minute in range(0, lt_sec + 60, 60):
                    dt_minute = datetime.time(hour=int(minute/3600), minute=int((minute % 3600)/60))
                    staytime_dist[dt_minute] += 1

                for minute in range(at_sec, 86400, 60):
                    dt_minute = datetime.time(hour=int(minute/3600), minute=int((minute % 3600)/60))
                    staytime_dist[dt_minute] += 1

        if arriving_leaving == 'f':
            tid = points_in_loc[pid1][4]

            next_loc = traj_from_to_loc[tid][1]

            nextloc_count[next_loc] += 1

            ts1 = datetime.datetime.fromtimestamp(trajectories[tid].start_point()[2])
            ts2 = datetime.datetime.fromtimestamp(trajectories[tid].end_point()[2])

            at = ts1.replace(second=0, microsecond=0)
            lt = ts2.replace(second=0, microsecond=0)

            midnight_at = at.replace(hour=0, minute=0)
            midnight_lt = lt.replace(hour=0, minute=0)

            at_sec = int((at - midnight_at).total_seconds())
            lt_sec = int((lt - midnight_lt).total_seconds())

            for minute in range(at_sec, lt_sec + 60, 60):
                dt_minute = datetime.time(hour=int(minute/3600), minute=int((minute % 3600)/60))
                nextloc_dist[next_loc][dt_minute] += 1

    spatial_points = list()
    for p in list(points_in_loc.values()):
        spatial_points.append(p[0:2])

    loc_rg = radius_of_gyration(spatial_points, location_prototype, spherical_distance)
    loc_entropy = entropy(list(nextloc_count.values()), classes=len(nextloc_count))

    res = {
        'staytime_dist': default_to_regular(staytime_dist),
        'nextloc_dist': default_to_regular(nextloc_dist),
        'nextloc_count': nextloc_count,
        'loc_rg': loc_rg,
        'loc_entropy': loc_entropy,
    }

    return res


def get_locations_features(points, traj_from_to_loc, location_points, location_prototype, imh):
    res = dict()

    for lid in location_points:

        points_in_loc = dict()
        for pid in location_points[lid]:
            points_in_loc[pid] = points[pid]

        lf_res = get_location_features(points_in_loc, traj_from_to_loc, location_prototype[lid], imh)
        lf_res['loc_support'] = len(location_points[lid])
        res[lid] = lf_res

    staytime_tot_dist = defaultdict(list)
    for lid in res:
        staytime_dist = res[lid]['staytime_dist']

        for ts in staytime_dist:
            staytime_tot_dist[ts].append(staytime_dist[ts])

    staytime_totals = dict()
    for ts in staytime_tot_dist:
        staytime_totals[ts] = np.sum(staytime_tot_dist[ts])

    for lid in res:
        staytime_dist = res[lid]['staytime_dist']
        staytime_ndist = dict()
        for ts in staytime_dist:
            staytime_ndist[ts] = 1.0 * staytime_dist[ts] / staytime_totals[ts]
        res[lid]['staytime_ndist'] = staytime_ndist

    return res


# Specific code for VFI pilot

def get_locations_durations_VFI(trajectories, traj_from_to_loc, loc_ft_mid):

    loc_arrive_leave_list = defaultdict(list)
    # last_loc_leave = None
    last_loc_arrive = None
    # last_time_leave = None
    last_time_arrive = None

    loc_timeslot_set_andrienko = dict()
    mov_timeslot_set_andrienko = dict()

    sorted_traj = sorted(trajectories, key=lambda t: trajectories[t].start_point()[2])

    for tid in sorted_traj:
        traj = trajectories[tid]
        loc_leave = traj_from_to_loc[tid][0]
        loc_arrive = traj_from_to_loc[tid][1]
        time_leave = datetime.datetime.fromtimestamp(traj.start_point()[2])
        time_arrive = datetime.datetime.fromtimestamp(traj.end_point()[2])

        if last_loc_arrive is not None:
            if loc_leave == last_loc_arrive:
                loc_arrive_leave_list[last_loc_arrive].append([last_time_arrive, time_leave])

            for ts in pd.date_range(last_time_arrive, time_leave, freq='H'):
                if loc_leave not in loc_timeslot_set_andrienko:
                    loc_timeslot_set_andrienko[last_loc_arrive] = {'weekdays': defaultdict(set), 'weekend': defaultdict(set)}
                type_of_day = 'weekend' if ts.weekday() >= 5 else 'weekdays'
                timeslot = get_timeslot(ts)
                day_key = (ts.year, ts.month, ts.day)
                loc_timeslot_set_andrienko[last_loc_arrive][type_of_day][timeslot].add(day_key)

        last_loc_arrive = loc_arrive
        last_time_arrive = time_arrive
        # print(tid, loc_leave, loc_arrive, time_leave, time_arrive)

        if loc_leave not in loc_timeslot_set_andrienko:
            loc_timeslot_set_andrienko[loc_leave] = {'weekdays': defaultdict(set), 'weekend': defaultdict(set)}
        type_of_day = 'weekend' if time_leave.weekday() >= 5 else 'weekdays'
        timeslot = get_timeslot(time_leave)
        day_key = (time_leave.year, time_leave.month, time_leave.day)
        loc_timeslot_set_andrienko[loc_leave][type_of_day][timeslot].add(day_key)

        if loc_arrive not in loc_timeslot_set_andrienko:
            loc_timeslot_set_andrienko[loc_arrive] = {'weekdays': defaultdict(set), 'weekend': defaultdict(set)}
        type_of_day = 'weekend' if time_arrive.weekday() >= 5 else 'weekdays'
        timeslot = get_timeslot(time_arrive)
        day_key = (time_arrive.year, time_arrive.month, time_arrive.day)
        loc_timeslot_set_andrienko[loc_arrive][type_of_day][timeslot].add(day_key)

        mid = loc_ft_mid[(loc_leave, loc_arrive)]
        if mid not in mov_timeslot_set_andrienko:
            mov_timeslot_set_andrienko[mid] = {'weekdays': defaultdict(set), 'weekend': defaultdict(set)}
        type_of_day = 'weekend' if time_leave.weekday() >= 5 else 'weekdays'
        timeslot = get_timeslot(time_leave)
        day_key = (time_leave.year, time_leave.month, time_leave.day)
        mov_timeslot_set_andrienko[mid][type_of_day][timeslot].add(day_key)

    stays_at_locations = []
    loc_arrive_leave_list = default_to_regular(loc_arrive_leave_list)
    for loc, arrive_leave_list in loc_arrive_leave_list.items():
        for al in arrive_leave_list:
            print("Loc dur ", loc, al[1], al[0], (al[1] - al[0]).total_seconds())
            stay_dict = {'location': str(loc), 'arrival_time': al[0], 'departure_time': al[1],
                         'stay_duration': (al[1] - al[0]).total_seconds()}
            stays_at_locations.append(stay_dict)

    print("stays")
    print(stays_at_locations)
    # for k, v in staytime_durations.items():
    #     print(k, v)

    # res = dict()
    # for loc, stay_info in staytime_durations.items():
    #     res[loc] = {
    #         'arrival_time': stay_info
    #         'avg_staytime': np.mean(stay_info[2]),
    #     }
    #
    return stays_at_locations


# End Specific code for VFI


def interquartile_filter(x):
    q1 = np.percentile(x, 25)
    q3 = np.percentile(x, 75)
    iqr = q3 - q1
    y = list()
    for x0 in x:
        if q1 - 1.5 * iqr <= x0 <= q3 + 1.5 * iqr:
            y.append(x0)
    return y


def get_movements_features(movement_traj, imh):

    trajectories = imh['trajectories']

    res = dict()
    for mid in movement_traj:
        traj_in_movement = movement_traj[mid][1]
        movement_support = len(traj_in_movement)
        movement_lengths = list()
        movement_durations = list()
        for tid in traj_in_movement:
            movement_lengths.append(trajectories[tid].length())
            movement_durations.append(trajectories[tid].duration())

        movement_lengths = interquartile_filter(movement_lengths)
        movement_durations = interquartile_filter(movement_durations)

        res[mid] = {
            'mov_support': movement_support,
            'typical_mov_length': np.median(movement_lengths),
            'avg_mov_length': np.mean(movement_lengths),
            'std_mov_length': np.std(movement_lengths),
            'typical_mov_duration': datetime.timedelta(seconds=np.median(movement_durations)),
            'avg_mov_duration': datetime.timedelta(seconds=np.mean(movement_durations)),
            'std_mov_duration': datetime.timedelta(seconds=np.std(movement_durations)),
        }

    return res


def get_movements_stats(movement_traj, regular_locs, imh):
    trajectories = imh['trajectories']

    movement_lengths = list()
    movement_durations = list()
    reg_movement_lengths = list()
    reg_movement_durations = list()
    reg_movs = dict()
    n_reg_traj = 0

    for mid in movement_traj:
        lft = movement_traj[mid][0]
        traj_in_movement = movement_traj[mid][1]
        for tid in traj_in_movement:
            movement_lengths.append(trajectories[tid].length())
            movement_durations.append(trajectories[tid].duration())
            if regular_locs is not None:
                if lft[0] in regular_locs and lft[1] in regular_locs:
                    reg_movs[mid] = 0
                    reg_movement_lengths.append(trajectories[tid].length())
                    reg_movement_durations.append(trajectories[tid].duration())
                    n_reg_traj += 1

    movement_lengths = interquartile_filter(movement_lengths)
    movement_durations = interquartile_filter(movement_durations)

    if len(reg_movement_lengths) > 0:
        reg_movement_lengths = interquartile_filter(reg_movement_lengths)
        reg_movement_durations = interquartile_filter(reg_movement_durations)

    avg_mov_duration = datetime.timedelta(seconds=np.mean(movement_durations))
    std_mov_duration = datetime.timedelta(seconds=np.std(movement_durations))

    if len(reg_movement_lengths) > 0:
        avg_reg_mov_duration = datetime.timedelta(seconds=np.mean(reg_movement_durations))
        std_reg_mov_duration = datetime.timedelta(seconds=np.std(reg_movement_durations))
    else:
        avg_reg_mov_duration = avg_mov_duration
        std_reg_mov_duration = std_mov_duration

    res = {
        'n_reg_movs': len(reg_movs),
        'avg_mov_length': np.mean(movement_lengths),
        'std_mov_length': np.std(movement_lengths),
        'avg_mov_duration': avg_mov_duration,
        'std_mov_duration': std_mov_duration,
        'avg_reg_mov_length': np.mean(reg_movement_lengths),
        'std_reg_mov_length': np.std(reg_movement_lengths),
        'avg_reg_mov_duration': avg_reg_mov_duration,
        'std_reg_mov_duration': std_reg_mov_duration,
        'n_reg_traj': n_reg_traj,
    }

    return res


def closest_point_on_segment_minsup(a, b, p):
    sx1 = a[0]
    sx2 = b[0]
    sy1 = a[1]
    sy2 = b[1]
    px = p[0]
    py = p[1]

    x_delta = sx2 - sx1
    y_delta = sy2 - sy1

    if x_delta == 0 and y_delta == 0:
        return p

    u = ((px - sx1) * x_delta + (py - sy1) * y_delta) / (x_delta * x_delta + y_delta * y_delta)

    if u < 0.00001:
        closest_point = a
    elif u > 1:
        closest_point = b
    else:
        cp_x = sx1 + u * x_delta
        cp_y = sy1 + u * y_delta
        closest_point = [cp_x, cp_y]

    return closest_point


def get_minimum_support(locations_support):
    x = []
    y = []

    sorted_support = sorted(locations_support)

    for i, s in enumerate(sorted_support):
        x.append(1.0 * i)
        y.append(1.0 * s)

    # plt.plot(x, y)
    # plt.plot([0, len(x) - 1], [y[0], y[len(y)-1]])

    max_d = -float('infinity')
    index = 0

    a = [x[0], y[0]]
    b = [x[len(x)-1], y[len(y)-1]]

    for i in range(0, len(x)):
        p = [x[i], y[i]]
        c = closest_point_on_segment_minsup(a, b, p)
        d = math.sqrt((c[0]-x[i])**2 + (c[1]-y[i])**2)

        if d > max_d:
            max_d = d
            index = i

    return sorted_support[index]


def detect_regular_locations(location_points, loc_nextlocs):

    loc_support = dict()

    for lid in location_points:
        loc_support[lid] = len(location_points[lid])

    loc_min_sup = get_minimum_support(list(loc_support.values()))

    regular_locs = dict()
    for lid in loc_support:
        if loc_support[lid] >= loc_min_sup:
            regular_locs[lid] = loc_nextlocs[lid]

    is_dag = False
    while not is_dag and not len(regular_locs) <= 2:
        is_dag = True
        for lid in regular_locs:
            has_an_out_mov = False
            for lid_out in regular_locs[lid]:
                if lid_out in regular_locs and lid_out != lid:
                    has_an_out_mov = True
                    break
            if not has_an_out_mov:
                del regular_locs[lid]
                is_dag = False
                break

    return regular_locs, loc_min_sup, is_dag


def caclulate_regular_rgen(regular_locs, loc_res, points):

    spatial_points = list()
    rloc_support = dict()
    for rlid in regular_locs:
        rloc_support[rlid] = len(loc_res['location_points'][rlid])
        for pid in loc_res['location_points'][rlid]:
            p = points[pid]
            spatial_points.append(p[0:2])

    cm = np.mean(spatial_points, axis=0)
    rrg = radius_of_gyration(spatial_points, cm, spherical_distance)
    ren = entropy(list(rloc_support.values()), classes=len(rloc_support))

    return rrg, ren


def get_events_features(events, traj_location_from_to):

    loc_from_to_info = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    loc_from_to_features = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

    for tid in traj_location_from_to:
        if tid not in events:
            continue

        lf, lt = traj_location_from_to[tid]
        lft = (lf, lt)
        for event in events[tid]:
            evnt_type = event['event_type']
            speed = event['speed']
            max_acc = event['max_acc']
            avg_acc = event['avg_acc']
            duration = event['duration']
            # location_type = event['location_type']

            loc_from_to_info[lft][evnt_type]['speed'].append(speed)
            loc_from_to_info[lft][evnt_type]['max_acc'].append(max_acc)
            loc_from_to_info[lft][evnt_type]['avg_acc'].append(avg_acc)
            loc_from_to_info[lft][evnt_type]['duration'].append(duration)

    for lft in loc_from_to_info:
        for evnt_type in loc_from_to_info[lft]:
            speed = loc_from_to_info[lft][evnt_type]['speed']
            loc_from_to_features[lft][evnt_type]['speed'] = [np.mean(speed), np.std(speed)]
            max_acc = loc_from_to_info[lft][evnt_type]['max_acc']
            loc_from_to_features[lft][evnt_type]['max_acc'] = [np.mean(max_acc), np.std(max_acc)]
            avg_acc = loc_from_to_info[lft][evnt_type]['avg_acc']
            loc_from_to_features[lft][evnt_type]['avg_acc'] = [np.mean(avg_acc), np.std(avg_acc)]
            duration = loc_from_to_info[lft][evnt_type]['duration']
            loc_from_to_features[lft][evnt_type]['duration'] = [np.mean(duration), np.std(duration)]
            loc_from_to_features[lft][evnt_type]['count'] = [len(speed)]

    loc_from_to_features = default_to_regular(loc_from_to_features)

    return loc_from_to_features


def get_timeslot(time_loc):
    if datetime.time(6) < time_loc.time() <= datetime.time(9):
        return '06-09'
    elif datetime.time(9) < time_loc.time() <= datetime.time(18):
        return '09-18'
    elif datetime.time(18) < time_loc.time() <= datetime.time(21):
        return '18-21'
    else:
        return '21-06'


def get_locations_durations(trajectories, traj_from_to_loc, loc_ft_mid):

    loc_arrive_leave_list = defaultdict(list)
    # last_loc_leave = None
    last_loc_arrive = None
    # last_time_leave = None
    last_time_arrive = None

    loc_timeslot_set_andrienko = dict()
    mov_timeslot_set_andrienko = dict()

    sorted_traj = sorted(trajectories, key=lambda t: trajectories[t].start_point()[2])

    for tid in sorted_traj:
        traj = trajectories[tid]
        loc_leave = traj_from_to_loc[tid][0]
        loc_arrive = traj_from_to_loc[tid][1]
        time_leave = datetime.datetime.fromtimestamp(traj.start_point()[2])
        time_arrive = datetime.datetime.fromtimestamp(traj.end_point()[2])

        if last_loc_arrive is not None:
            if loc_leave == last_loc_arrive:
                loc_arrive_leave_list[last_loc_arrive].append([last_time_arrive, time_leave])

            for ts in pd.date_range(last_time_arrive, time_leave, freq='H'):
                if loc_leave not in loc_timeslot_set_andrienko:
                    loc_timeslot_set_andrienko[last_loc_arrive] = {'weekdays': defaultdict(set), 'weekend': defaultdict(set)}
                type_of_day = 'weekend' if ts.weekday() >= 5 else 'weekdays'
                timeslot = get_timeslot(ts)
                day_key = (ts.year, ts.month, ts.day)
                loc_timeslot_set_andrienko[last_loc_arrive][type_of_day][timeslot].add(day_key)

        last_loc_arrive = loc_arrive
        last_time_arrive = time_arrive
        # print(tid, loc_leave, loc_arrive, time_leave, time_arrive)

        if loc_leave not in loc_timeslot_set_andrienko:
            loc_timeslot_set_andrienko[loc_leave] = {'weekdays': defaultdict(set), 'weekend': defaultdict(set)}
        type_of_day = 'weekend' if time_leave.weekday() >= 5 else 'weekdays'
        timeslot = get_timeslot(time_leave)
        day_key = (time_leave.year, time_leave.month, time_leave.day)
        loc_timeslot_set_andrienko[loc_leave][type_of_day][timeslot].add(day_key)

        if loc_arrive not in loc_timeslot_set_andrienko:
            loc_timeslot_set_andrienko[loc_arrive] = {'weekdays': defaultdict(set), 'weekend': defaultdict(set)}
        type_of_day = 'weekend' if time_arrive.weekday() >= 5 else 'weekdays'
        timeslot = get_timeslot(time_arrive)
        day_key = (time_arrive.year, time_arrive.month, time_arrive.day)
        loc_timeslot_set_andrienko[loc_arrive][type_of_day][timeslot].add(day_key)

        mid = loc_ft_mid[(loc_leave, loc_arrive)]
        if mid not in mov_timeslot_set_andrienko:
            mov_timeslot_set_andrienko[mid] = {'weekdays': defaultdict(set), 'weekend': defaultdict(set)}
        type_of_day = 'weekend' if time_leave.weekday() >= 5 else 'weekdays'
        timeslot = get_timeslot(time_leave)
        day_key = (time_leave.year, time_leave.month, time_leave.day)
        mov_timeslot_set_andrienko[mid][type_of_day][timeslot].add(day_key)

    staytime_durations = defaultdict(list)
    loc_arrive_leave_list = default_to_regular(loc_arrive_leave_list)
    for loc, arrive_leave_list in loc_arrive_leave_list.items():
        for al in arrive_leave_list:
            staytime_durations[loc].append((al[1] - al[0]).total_seconds())

    staytime_durations = default_to_regular(staytime_durations)
    # for k, v in staytime_durations.items():
    #     print(k, v)

    res = dict()
    for loc, sd in staytime_durations.items():
        res[loc] = {
            #'staytime_durations_list': sd,
            'avg_staytime': np.mean(sd),
            'std_staytime': np.std(sd),
        }

    loc_timeslot_count_andrienko = dict()
    for loc in loc_timeslot_set_andrienko:
        loc_timeslot_count_andrienko[loc] = {'weekdays': dict(), 'weekend': dict()}
        for type_of_day in loc_timeslot_set_andrienko[loc]:
            for timeslot in loc_timeslot_set_andrienko[loc][type_of_day]:
                val = len(loc_timeslot_set_andrienko[loc][type_of_day][timeslot])
                loc_timeslot_count_andrienko[loc][type_of_day][timeslot] = val

    mov_timeslot_count_andrienko = dict()
    for mid in mov_timeslot_set_andrienko:
        mov_timeslot_count_andrienko[mid] = {'weekdays': dict(), 'weekend': dict()}
        for type_of_day in mov_timeslot_set_andrienko[mid]:
            for timeslot in mov_timeslot_set_andrienko[mid][type_of_day]:
                val = len(mov_timeslot_set_andrienko[mid][type_of_day][timeslot])
                mov_timeslot_count_andrienko[mid][type_of_day][timeslot] = val

    return res, loc_timeslot_count_andrienko, mov_timeslot_count_andrienko


def build_imn(imh, reg_loc=True, events=None, verbose=False, pilot="SIS"):

    n_traj = len(imh['trajectories'])
    if verbose:
        print(datetime.datetime.now(), 'user %s - trajectories %s' % (imh['uid'], n_traj))

    if verbose:
        print(datetime.datetime.now(), 'extract start end')
    points, traj_from_to = get_points_trajfromto(imh)

    # added traj s/e times
    tid_se_times = {}
    for tid in traj_from_to:
        start_time = int(imh['trajectories'][tid].start_point()[2])
        end_time = int(imh['trajectories'][tid].end_point()[2])
        tid_se_times[tid] = (start_time, end_time)

    start_time = time.time()
    if verbose:
        print(datetime.datetime.now(), 'location detection')
    loc_res = locations_detection(points)
    print("--- %s seconds Location Extraction ---" % (time.time() - start_time))

    if loc_res is None:
        return None

    if verbose:
        print(datetime.datetime.now(), 'movement detection')
    mov_res = movements_detection(loc_res['pid_lid'], traj_from_to, imh)
    # print(traj_from_to)
    # print(mov_res)

    n_locs = len(loc_res['location_points'])
    n_movs = len(mov_res['movement_traj'])

    rrg = ren = loc_min_sup = None
    regular_locs = loc_res['location_points']
    n_reg_locs = len(regular_locs)

    if reg_loc:
        if verbose:
            print(datetime.datetime.now(), 'regularities detection')
        regular_locs, loc_min_sup, is_dag = detect_regular_locations(
            loc_res['location_points'], mov_res['loc_nextlocs'])
        n_reg_locs = len(regular_locs)
        rrg, ren = caclulate_regular_rgen(regular_locs, loc_res, points)

    if verbose:
        print(datetime.datetime.now(), 'location features extraction')


    lf_res = get_locations_features(points, mov_res['traj_from_to_loc'], loc_res['location_points'],
                                    loc_res['location_prototype'], imh)

    lf_dur, lf_ac, mf_ac = get_locations_durations(imh['trajectories'], mov_res['traj_from_to_loc'], mov_res['lft_mid'])

    for k, v in lf_res.items():
        if k in lf_dur:
            lf_res[k].update(lf_dur[k])
    for k, v in lf_res.items():
        if k in lf_ac:
            lf_res[k]['timeslot_count'] = lf_ac[k]

    if verbose:
        print(datetime.datetime.now(), 'movement features extraction')
    mf_res = get_movements_features(mov_res['movement_traj'], imh)
    ms_res = get_movements_stats(mov_res['movement_traj'], regular_locs, imh)
    for k, v in mf_res.items():
        if k in mf_ac:
            mf_res[k]['timeslot_count'] = mf_ac[k]

    if pilot == "SIS":
        imn = {
            'point_location': loc_res['pid_lid'],
            'traj_points_from_to': traj_from_to,
            'traj_location_from_to': mov_res['traj_from_to_loc'],
            'tid_se_times': tid_se_times,

            'location_points': default_to_regular(loc_res['location_points']),
            'regular_locations': list(regular_locs.keys()),
            'location_prototype': loc_res['location_prototype'],
            'location_nextlocs': mov_res['loc_nextlocs'],
            'location_features': lf_res,

            'movement_traj': mov_res['movement_traj'],
            'movement_prototype': mov_res['movement_prototype'],
            'location_from_to_movement': mov_res['lft_mid'],
            'mov_features': mf_res,

            'n_traj': n_traj,
            'n_reg_traj': ms_res['n_reg_traj'],
            'n_locs': n_locs,
            'n_reg_locs': n_reg_locs,
            'n_movs': n_movs,
            'n_reg_movs': ms_res['n_reg_movs'],
            'rg': loc_res['rg'],
            'rrg': rrg,
            'entropy': loc_res['entropy'],
            'rentropy': ren,
            'avg_mov_length': ms_res['avg_mov_length'],
            'std_mov_length': ms_res['std_mov_length'],
            'avg_mov_duration': ms_res['avg_mov_duration'],
            'std_mov_duration': ms_res['std_mov_duration'],
            'avg_reg_mov_length': ms_res['avg_reg_mov_length'],
            'std_reg_mov_length': ms_res['std_reg_mov_length'],
            'avg_reg_mov_duration': ms_res['avg_reg_mov_duration'],
            'std_reg_mov_duration': ms_res['std_reg_mov_duration'],
            'loc_tosca_cut': loc_res['loc_tosca_cut'],
            'loc_sup_cut': loc_min_sup,

            'graph': mov_res['graph'],
        }
    if pilot == "VFI":
        stays_VFI = get_locations_durations_VFI(imh['trajectories'], mov_res['traj_from_to_loc'], mov_res['lft_mid'])
        imn = {
            'locations': loc_res['location_prototype'],
            'stays': stays_VFI,
            #'location_features': lf_res,
            #'location_nextlocs': mov_res['loc_nextlocs'],
            'location_from_to_movement': mov_res['lft_mid'],
            'mov_features': mf_res,

            'graph': mov_res['graph'],
        }

    if events is not None:
        traj_location_from_to = imn['traj_location_from_to']
        evnt_res = get_events_features(events, traj_location_from_to)
        imn['events'] = evnt_res

    return imn


def get_lonlat(location):
    lon = location[0]
    lat = location[1]
    return lon, lat
