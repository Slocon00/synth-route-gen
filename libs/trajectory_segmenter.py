import math
import numpy as np

from scipy import stats

from trajectory import Trajectory
from mobility_distance_functions import spherical_distance


def thompson_test(data, point, alpha=0.05):
    sd = np.std(data)
    mean = np.mean(data)
    delta = abs(point - mean)

    t_a2 = stats.t.ppf(1-(alpha/2.), len(data)-2)
    tau = (t_a2 * (len(data)-1)) / (math.sqrt(len(data)) * math.sqrt(len(data)-2 + t_a2**2))

    return delta > tau * 2*sd


def segment_trajectories(alltraj, uid, temporal_thr=60, spatial_thr=50, max_speed=0.07):
    # temporal_thr = 120 # seconds
    # spatial_thr = 50 # meters
    # max_speed = 0.07 # km/s
    max_speed = max_speed * 1000
    spatial_thr = spatial_thr

    #traj_list = list()
    traj_list = dict()

    tid = 0
    traj = list()
    is_a_new_traj = True
    p = None
    length = 0.0
    ref_p = None  # for stop detection
    first_iteration = True

    p_index = 0
    next_p_index = 0
    ref_p_index = 0

    for i in range(0, len(alltraj)):

        next_p = alltraj[i]
        next_p_index = i

        if first_iteration:  # first iteration
            p = next_p
            p_index = next_p_index
            ref_p = p  # for stop detection
            traj = [p]
            length = 0.0
            is_a_new_traj = True
            first_iteration = False
        else:                 # all the others
            spatial_dist = spherical_distance(p, next_p)
            temporal_dist = next_p[2] - p[2]

            ref_distance = spherical_distance(ref_p, next_p)  # for stop detection
            ref_time_diff = next_p[2] - ref_p[2]  # for stop detection

            # Ignore extreme jump (with speed > 250km/h = 0.07km/s)
            if ref_distance > max_speed * ref_time_diff:
                continue

            if temporal_dist > temporal_thr or (ref_time_diff > temporal_thr and ref_distance < spatial_thr):
                # ended trajectory (includes case with long time gap)
                if len(traj) > 1 and not is_a_new_traj:
                    start_time = traj[0][2]
                    end_time = traj[-1][2]
                    duration = end_time - start_time
                    trajToInsert = Trajectory(id=tid, object=traj, vehicle=uid,
                                              length=length, duration=duration,
                                              start_time=start_time, end_time=end_time)
                    traj_list[tid] = trajToInsert
                    #traj_list.append()

                # Create a new trajectory
                traj = [p[:2] + [next_p[2]]]  # 1st fake point with last position previous traj and new timestamp
                p = next_p
                p_index = next_p_index
                ref_p = p  # for stop detection
                ref_p_index = p_index
                traj.append(p)
                length = 0.0
                tid += 1
                is_a_new_traj = True
            else:
                if is_a_new_traj and len(traj) > 1:
                    traj[1] = [traj[1][0], traj[1][1], traj[0][2] +int((next_p[2] - traj[0][2])/2)]
                is_a_new_traj = False

                p = next_p
                p_index = next_p_index
                if ref_distance > spatial_thr:
                    ref_p = p  # for stop detection
                    ref_p_index = p_index
                traj.append(p)
                length += spatial_dist

    if len(traj) > 1 and not is_a_new_traj:
        start_time = traj[0][2]
        end_time = traj[-1][2]
        duration = end_time - start_time
        trajToInsert = Trajectory(id=tid, object=traj, vehicle=uid, length=length, duration=duration,
                                  start_time=start_time, end_time=end_time)
        traj_list[tid] = trajToInsert
        #traj_list.append()

    return traj_list


def segment_trajectories_random(alltraj, uid, nbr_traj_min=None, nbr_traj_max=None, nbr_traj=None):
    nbr_traj_min = 2 if nbr_traj_min is None else nbr_traj_min
    nbr_traj_max = int(len(alltraj) / 2) if nbr_traj_max is None else nbr_traj_max

    nbr_traj = np.random.randint(nbr_traj_min, nbr_traj_max + 1) if nbr_traj is None else nbr_traj
    new_traj_marker = int(len(alltraj) / nbr_traj)

    traj_list = list()

    tid = 0
    traj = list()
    is_a_new_traj = True
    p = None
    first_iteration = True
    length = 0.0

    for i in range(0, len(alltraj)):

        next_p = alltraj[i]

        if first_iteration:  # first iteration
            p = next_p
            traj = [p]
            length = 0.0
            is_a_new_traj = True
            first_iteration = False
        else:  # all the others
            spatial_dist = spherical_distance(p, next_p)
            if i % new_traj_marker == 0 and new_traj_marker > 0:
                if len(traj) > 1 and not is_a_new_traj:
                    start_time = traj[0][2]
                    end_time = traj[-1][2]
                    duration = end_time - start_time
                    traj_list.append(Trajectory(id=tid, object=traj, vehicle=uid,
                                                length=length, duration=duration,
                                                start_time=start_time, end_time=end_time))

                # Create a new trajectory
                p = next_p
                traj = [p]
                length = 0.0
                tid += 1
                is_a_new_traj = True
            else:
                is_a_new_traj = False
                p = next_p
                traj.append(p)
                length += spatial_dist

    if len(traj) > 1 and not is_a_new_traj:
        start_time = traj[0][2]
        end_time = traj[-1][2]
        duration = end_time - start_time
        traj_list.append(Trajectory(id=tid, object=traj, vehicle=uid, length=length, duration=duration,
                                    start_time=start_time, end_time=end_time))

    return traj_list


def segment_trajectories_random2(alltraj, uid, nbr_traj):
    new_traj_marker = int(len(alltraj) / nbr_traj)

    traj_list = list()

    tid = 0
    traj = list()
    is_a_new_traj = True
    p = None
    first_iteration = True
    length = 0.0

    for i in range(0, len(alltraj)):

        next_p = alltraj[i]

        if first_iteration:  # first iteration
            p = next_p
            traj = [p]
            length = 0.0
            is_a_new_traj = True
            first_iteration = False
        else:  # all the others
            spatial_dist = spherical_distance(p, next_p)
            if i % new_traj_marker == 0 and new_traj_marker > 0:
                if len(traj) > 1 and not is_a_new_traj:
                    start_time = traj[0][2]
                    end_time = traj[-1][2]
                    duration = end_time - start_time
                    traj_list.append(Trajectory(id=tid, object=traj, vehicle=uid,
                                                length=length, duration=duration,
                                                start_time=start_time, end_time=end_time))

                # Create a new trajectory
                p = next_p
                traj = [p]
                length = 0.0
                tid += 1
                is_a_new_traj = True
            else:
                is_a_new_traj = False
                p = next_p
                traj.append(p)
                length += spatial_dist

    if len(traj) > 1 and not is_a_new_traj:
        start_time = traj[0][2]
        end_time = traj[-1][2]
        duration = end_time - start_time
        traj_list.append(Trajectory(id=tid, object=traj, vehicle=uid, length=length, duration=duration,
                                    start_time=start_time, end_time=end_time))

    return traj_list



def get_stop_times(traj_list):
    stop_time_list = list()
    for i in range(1, len(traj_list)):
        arrival_time = traj_list[i-1].end_time()
        leaving_time = traj_list[i].start_time()
        stop_time_list.append(leaving_time - arrival_time)
    return stop_time_list


def moving_avg(a, w=3):
    ma = list()
    for i in range(0, len(a) - w):
        m = np.mean(a[i:i+w])
        ma.append(m)
    return ma


def moving_median(a, w=3):
    ma = list()
    for i in range(0, len(a) - w):
        m = np.median(a[i:i+w])
        ma.append(m)
    return ma


def segment_trajectories_user_adaptive(alltraj, uid, temporal_thr=60, spatial_thr=50, max_speed=0.07,
                                       gap=60, max_lim=3600*48, window=3, smooth_fun=moving_avg, min_size=10,
                                       return_cut=False):
    traj_list = segment_trajectories(alltraj, uid, temporal_thr, spatial_thr, max_speed)
    stop_time_list = get_stop_times(traj_list)

    time_stop_values = range(gap, max_lim + gap, gap)
    # time_stop_values = np.concatenate([np.arange(60, 60*10, 60), np.arange(60*10, max_lim + gap, gap)])
    stop_time_count, _ = np.histogram(stop_time_list, bins=time_stop_values)

    stop_time_count_ma = smooth_fun(stop_time_count[::-1], window)
    time_stop_values_ma = smooth_fun(time_stop_values[::-1], window)

    user_temporal_thr = 1200
    for cut in range(len(stop_time_count_ma) - 1, min_size, -1):
        if thompson_test(stop_time_count_ma[:cut], stop_time_count_ma[cut]):
            #  and time_stop_values_ma[cut] > time_stop_values_ma[cut]:
            # print(cut, stop_time_count_ma[cut], time_stop_values_ma[cut])
            user_temporal_thr = time_stop_values_ma[cut]
            break

    traj_list = segment_trajectories(alltraj, uid, user_temporal_thr, spatial_thr, max_speed)
    if return_cut:
        return traj_list, user_temporal_thr
    return traj_list
