import json

from mobility_distance_functions import *


__author__ = 'Riccardo Guidotti'


class Trajectory:
    """description"""

    def __init__(self, id, object, vehicle, length=None, duration=None, start_time=None, end_time=None):
        self.id = id
        self.object = object
        self.vehicle = vehicle
        self._length = length
        self._duration = duration
        self._start_time = start_time
        self._end_time = end_time

    def id(self):
        return self.id

    def object(self):
        return self.object

    def vehicle(self):
        return self.vehicle

    def num_points(self):
        return self.__len__()

    def point_n(self, n):
        return self.object[n]

    def start_point(self):
        return self.point_n(0)

    def end_point(self):
        return self.point_n(len(self)-1)

    def start_time(self):
        return self._start_time

    def end_time(self):
        return self._end_time

    def length(self):
        if self._length is None:
            length = 0
            for i in range(0, len(self.object)-1, 1):
                p1 = self.point_n(i)
                p2 = self.point_n(i+1)
                dist = spherical_distance(p1, p2)
                length += dist
            self._length = length
        return self._length

    def longest_segment(self):
        # restituisce la distanza maggiore tra due punti consecutivi della traiettoria
        longest_dist = 0
        for i in range(0, len(self.object)-1):
            p1 = self.point_n(i)
            p2 = self.point_n(i+1)
            dist = spherical_distance(p1, p2)
            if dist > longest_dist:
                longest_dist = dist
        return longest_dist

    def duration(self):
        if self._duration is None:
            duration = 0
            for i in range(0, len(self.object) - 1, 1):
                p1 = self.point_n(i)
                p2 = self.point_n(i + 1)
                dist = p2[2] - p1[2]
                duration += dist
            self._duration = duration
        return self._duration

    def __len__(self):
        return len(self.object)

    def __str__(self):
        return json.dumps({'id': self.id, 'vehicle': self.vehicle, 'object': self.object})

    def to_json(self):
        return self.__dict__

    # def __repr__(self):
    #     return self.to_json()

    # def __unicode__(self):
    #     return self.__str__()



def calculate_traj_approximation(traj1, traj2, pred_thr, last_prop, time_mod=86400):
    """
    Calculate the approximation between the traejctory and the routine.
    """
    res = {
        'id_t1': traj1.id,
        'id_t2': traj2.id,
        'head': None,
        'tail': None,
        'dist': float('infinity')
    }

    # lool for the closest point on routine to the last point of traj
    t_last = traj1.end_point()
    min_dist = float('infinity')
    id_min = None
    for i in range(0, len(traj2) - 1):
        p = traj2.point_n(i)
        dist = spherical_distance(p, t_last)
        if dist < min_dist:
            min_dist = dist
            id_min = i

    cp = closest_point_on_segment(traj2.point_n(id_min), traj2.point_n(id_min + 1), t_last)

    # calcualte the distance between the two closest points
    dist = spherical_distance(cp, t_last)
    if last_prop == 0.0 and dist >= pred_thr:
        return res

    # cut the routine temporally from the beginning of traj to the time of the closest point
    t2 = cp[2] % time_mod
    t1 = traj2.start_point()[2] % time_mod
    traj2_cut = get_sub_trajectory(traj2, t1, t2)
    if traj2_cut is None or len(traj2_cut) < 3:
        return res

    # if the trajectory is shorter than the routine cut remove the initial part of the routine_cut
    if traj1.length() < traj2_cut.length():
        traj2_cut = get_sub_trajectory_keep_end(traj2_cut, traj1.length())

    # calculate the tail
    traj2_head = traj2_cut

    traj2_tail = get_sub_trajectory(traj2, t2, traj2.end_point()[2] % time_mod)

    if traj2_tail is None:
        traj2_tail = Trajectory(id=traj2.id, object=[traj2.end_point()], vehicle=traj2.vehicle)

    traj2_tail.object.insert(0, traj2_cut.end_point())

    if last_prop > 0.0:
        last_traj = get_sub_trajectory_keep_end(traj1, traj1.length() * last_prop)
        last_routine_head = get_sub_trajectory_keep_end(traj2_head, traj2_head.length() * last_prop)
        if len(last_traj) >= 2 and len(last_routine_head) >= 2:
            dist = trajectory_distance(last_traj, last_routine_head)
            if dist >= pred_thr:
                return res

    res['head'] = traj2_head
    res['tail'] = traj2_tail
    res['dist'] = dist

    return res


def get_sub_trajectory(traj, from_ts, to_ts, time_mod=86400):
    """
    Cut traj according to the temporal thresholds from and to.
    """
    t_start = traj.start_point()[2] % time_mod
    t_end = traj.end_point()[2] % time_mod

    if to_ts < t_start:
        return None

    if from_ts > t_end:
        return None

    if from_ts < t_start:
        from_ts = t_start

    if to_ts > t_end:
        to_ts = t_end

    id_sub = traj.id
    object_sub = list()
    vehicle_sub = traj.vehicle

    for i in range(0, len(traj)):
        ts = traj.point_n(i)[2] % time_mod
        if from_ts <= ts <= to_ts:
            object_sub.append(traj.point_n(i))
        if ts > to_ts:
            break

    sub_trajectory = Trajectory(id=id_sub, object=object_sub, vehicle=vehicle_sub)
    return sub_trajectory


def get_sub_trajectory_keep_end(traj, length):
    """
    Cut initial part of traj such that the length is respected.
    """
    if length >= traj.length():
        return traj

    id_sub = traj.id
    object_sub = [traj.end_point()]
    vehicle_sub = traj.vehicle

    tmp_length = 0
    for i in range(len(traj) - 1, 0, -1):
        p = traj.point_n(i)
        q = object_sub[len(object_sub) - 1]
        tmp_length += spherical_distance(q, p)
        if tmp_length >= length:
            break
        object_sub.insert(0, p)

    sub_trajectory = Trajectory(id=id_sub, object=object_sub, vehicle=vehicle_sub)
    return sub_trajectory


