import math
from geopy import distance

__author__ = 'Riccardo Guidotti'

def my_distance(coord1, coord2):
    new1 = [coord1[1], coord1[0]]
    new2 = [coord2[1], coord2[0]]
    dist = distance.distance(new1, new2).meters
    return dist

def my_normalization(data, minval, maxval):
    return (data-minval)/(maxval-minval)


def spherical_distance(a, b, measure='m'):
    lat1 = a[1]
    lon1 = a[0]
    lat2 = b[1]
    lon2 = b[0]
    if measure == 'km':
        R = 6367
    else:
        R = 6371000 #meters

    rlon1 = lon1 * math.pi / 180.0
    rlon2 = lon2 * math.pi / 180.0
    rlat1 = lat1 * math.pi / 180.0
    rlat2 = lat2 * math.pi / 180.0
    dlon = (rlon1 - rlon2) / 2.0
    dlat = (rlat1 - rlat2) / 2.0
    lat12 = (rlat1 + rlat2) / 2.0
    sindlat = math.sin(dlat)
    sindlon = math.sin(dlon)
    cosdlon = math.cos(dlon)
    coslat12 = math.cos(lat12)
    f = sindlat * sindlat * cosdlon * cosdlon + sindlon * sindlon * coslat12 * coslat12
    f = math.sqrt(f)
    f = math.asin(f) * 2.0 # the angle between the points
    f *= R
    return f


def start_end_distance(tr1, tr2):

    start1 = tr1.start_point()
    start2 = tr2.start_point()

    end1 = tr1.end_point()
    end2 = tr2.end_point()

    dist_start = spherical_distance(start1, start2)
    dist_end = spherical_distance(end1, end2)

    dist = dist_start + dist_end
    return dist





def point_at_time_agenda(a, b, ts):
    """
    Returns the points p at time ts between the points a and b
    """

    time_dist_a_b = b[2] - a[2]
    time_dist_a_p = ts

    # print time_dist_a_b, time_dist_a_p, '<<<<<'

    if time_dist_a_p >= time_dist_a_b:
        return b

    # find the distance from a to p
    # space_dist_a_b = spherical_distance(a, b)
    space_dist_a_b = math.sqrt((a[0] - b[0])**2 + (a[1] - b[1])**2)
    space_dist_a_p = 1.0 * time_dist_a_p / time_dist_a_b * space_dist_a_b

    # find point p
    p = [0, 0, a[2] + ts]

    if b[0] - a[0] == 0:
        return b

    m = (b[1] - a[1]) / (b[0] - a[0])
    p_0_1 = a[0] + space_dist_a_p / math.sqrt(1 + m**2)
    p_0_2 = a[0] - space_dist_a_p / math.sqrt(1 + m**2)
    p[0] = p_0_1 if p_0_1 > a[0] else p_0_2
    p[1] = m * (p[0] - a[0]) + a[1]

    return p


def point_at_time(a, b, ts, time_mod=86400):
    """
    Returns the points p at time ts between the points a and b
    """

    time_dist_a_b = (b[2]) % time_mod - (a[2]) % time_mod
    time_dist_a_p = ts

    # print (b[2] / 1000) % time_mod, (a[2] / 1000) % time_mod, ts
    # print time_dist_a_b, time_dist_a_p, '<<<<<'

    if time_dist_a_p >= time_dist_a_b:
        # print 'QUI'
        return b

    # find the distance from a to p
    # space_dist_a_b = spherical_distance(a, b)
    space_dist_a_b = math.sqrt((a[0] - b[0])**2 + (a[1] - b[1])**2)
    space_dist_a_p = 1.0 * time_dist_a_p / time_dist_a_b * space_dist_a_b

    # find point p
    p = [0, 0, a[2] + ts]

    if b[0] - a[0] == 0:
        # print 'QUO'
        return b

    m = (b[1] - a[1]) / (b[0] - a[0])
    p_0_1 = a[0] + space_dist_a_p / math.sqrt(1 + m**2)
    p_0_2 = a[0] - space_dist_a_p / math.sqrt(1 + m**2)
    p[0] = p_0_1 if p_0_1 > a[0] else p_0_2
    p[1] = m * (p[0] - a[0]) + a[1]

    return p


def __is_synch(p1, p2, time_th, time_mod=86400):
    ts1 = p1[2] % time_mod
    ts2 = p2[2] % time_mod
    return abs(ts1-ts2) >= time_th


def trajectory_distance_synch(tr1, tr2, time_th):
    if __is_synch(tr1.start_point(), tr2.start_point(), time_th) \
            and __is_synch(tr1.end_point(), tr2.end_point(), time_th):
        return float('infinity')
    else:
        return trajectory_distance(tr1, tr2)


def trajectory_distance_start_synch(tr1, tr2, time_th):
    if __is_synch(tr1.start_point(), tr2.start_point(), time_th):
        return float('infinity')
    else:
        return trajectory_distance(tr1, tr2)


def trajectory_distance_end_synch(tr1, tr2, time_th):
    if __is_synch(tr1.end_point(), tr2.end_point(), time_th):
        return float('infinity')
    else:
        return trajectory_distance(tr1, tr2)


def start_end_distance_synch(tr1, tr2, time_th):
    if __is_synch(tr1.start_point(), tr2.start_point(), time_th) \
            and __is_synch(tr1.end_point(), tr2.end_point(), time_th):
        return float('infinity')
    else:
        return start_end_distance(tr1, tr2)


def start_end_distance_start_synch(tr1, tr2, time_th):
    if __is_synch(tr1.start_point(), tr2.start_point(), time_th):
        return float('infinity')
    else:
        return start_end_distance(tr1, tr2)


def start_end_distance_end_synch(tr1, tr2, time_th):
    if __is_synch(tr1.end_point(), tr2.end_point(), time_th):
        return float('infinity')
    else:
        return start_end_distance(tr1, tr2)


def start_distance_synch(tr1, tr2, time_th):
    if __is_synch(tr1.start_point(), tr2.start_point(), time_th):
        return float('infinity')
    else:
        return start_distance(tr1, tr2)


def end_distance_synch(tr1, tr2, time_th):
    if __is_synch(tr1.end_point(), tr2.end_point(), time_th):
        return float('infinity')
    else:
        return end_distance(tr1, tr2)


def start_distance(tr1, tr2):

    start1 = tr1.start_point()
    start2 = tr2.start_point()

    dist_start = spherical_distance(start1, start2)

    dist = dist_start
    return dist


def end_distance(tr1, tr2):

    end1 = tr1.end_point()
    end2 = tr2.end_point()

    dist_end = spherical_distance(end1, end2)

    dist = dist_end
    return dist


def trajectory_distance(tr1, tr2):

    i1 = 0
    i2 = 0
    np = 0

    last_tr1 = tr1.point_n(i1)
    last_tr2 = tr2.point_n(i2)

    dist = spherical_distance(last_tr1, last_tr2)
    np += 1

    while True:

        step_tr1 = spherical_distance(last_tr1, tr1.point_n(i1+1))
        step_tr2 = spherical_distance(last_tr2, tr2.point_n(i2+1))

        if step_tr1 < step_tr2:
            i1 += 1
            last_tr1 = tr1.point_n(i1)
            last_tr2 = closest_point_on_segment(last_tr2, tr2.point_n(i2+1), last_tr1)
        elif step_tr1 > step_tr2:
            i2 += 1
            last_tr2 = tr2.point_n(i2)
            last_tr1 = closest_point_on_segment(last_tr1, tr1.point_n(i1+1), last_tr2)
        else:
            i1 += 1
            i2 += 1
            last_tr1 = tr1.point_n(i1)
            last_tr2 = tr2.point_n(i2)

        d = spherical_distance(last_tr1, last_tr2)

        dist += d
        np += 1

        if i1 >= (len(tr1)-1) or i2 >= (len(tr2)-1):
            break

    for i in range(i1, len(tr1)):
        d = spherical_distance(tr2.end_point(), tr1.point_n(i))
        dist += d
        np += 1

    for i in range(i2, len(tr2)):
        d = spherical_distance(tr1.end_point(), tr2.point_n(i))
        dist += d
        np += 1

    dist = 1.0 * dist / np

    return dist


def trajectory_distance2(tr1, tr2):

    i1 = 0
    i2 = 0
    np = 0

    last_tr1 = tr1.point_n(i1)
    last_tr2 = tr2.point_n(i2)

    tr1_length = tr1.length()
    tr2_length = tr2.length()

    tr_long = None
    tr_short = None
    if tr1_length <= tr2_length:
        """???"""

    dist = spherical_distance(last_tr1, last_tr2)
    np += 1

    while True:

        step_tr1 = spherical_distance(last_tr1, tr1.point_n(i1+1))
        step_tr2 = spherical_distance(last_tr2, tr2.point_n(i2+1))

        if step_tr1 < step_tr2:
            i1 += 1
            last_tr1 = tr1.point_n(i1)
            last_tr2 = closest_point_on_segment(last_tr2, tr2.point_n(i2+1), last_tr1)
        elif step_tr1 > step_tr2:
            i2 += 1
            last_tr2 = tr2.point_n(i2)
            last_tr1 = closest_point_on_segment(last_tr1, tr1.point_n(i1+1), last_tr2)
        else:
            i1 += 1
            i2 += 1
            last_tr1 = tr1.point_n(i1)
            last_tr2 = tr2.point_n(i2)

        d = spherical_distance(last_tr1, last_tr2)

        dist += d
        np += 1

        if i1 >= (len(tr1)-1) or i2 >= (len(tr2)-1):
            break

    for i in range(i1, len(tr1)):
        d = spherical_distance(tr2.end_point(), tr1.point_n(i))
        dist += d
        np += 1

    for i in range(i2, len(tr2)):
        d = spherical_distance(tr1.end_point(), tr2.point_n(i))
        dist += d
        np += 1

    dist = 1.0 * dist / np

    return dist


def closest_point_on_segment(a, b, p):
    sx1 = a[0]
    sx2 = b[0]
    sy1 = a[1]
    sy2 = b[1]
    sz1 = a[2]
    sz2 = b[2]
    px = p[0]
    py = p[1]

    x_delta = sx2 - sx1
    y_delta = sy2 - sy1
    z_delta = sz2 - sz1

    if x_delta == 0 and y_delta == 0:
        return p

    u = ((px - sx1) * x_delta + (py - sy1) * y_delta) / (x_delta * x_delta + y_delta * y_delta)
    if u < 0:
        closest_point = a
    elif u > 1:
        closest_point = b
    else:
        cp_x = sx1 + u * x_delta
        cp_y = sy1 + u * y_delta
        dist_a_cp = spherical_distance(a, [cp_x, cp_y, 0])
        if dist_a_cp != 0:
            cp_z = sz1 + z_delta / (spherical_distance(a, b) / spherical_distance(a, [cp_x, cp_y, 0]))
        else:
            cp_z = a[2]
        closest_point = [cp_x, cp_y, cp_z]

    return closest_point


def projection_factor(ls_s, ls_e, p):
    if p == ls_s:
        return 0.0
    if p == ls_e:
        return 1.0

    dx = ls_e[0] - ls_s[0]
    dy = ls_e[1] - ls_s[1]
    ls_len = dx * dx + dy * dy

    # handle zero-length segments
    if ls_len <= 0.0:
        return float("NaN")

    r = ( (p[0] - ls_s[0]) * dx + (p[1] - ls_s[1]) * dy ) / ls_len
    return r

def closest_point_on_segment2(ls_s, ls_e, p):
    factor = projection_factor(ls_s, ls_e, p)
    #print(factor)
    if (factor > 0 and factor < 1):
        if p == ls_s or p == ls_e:
            return (p, 0.0)
        cpx = ls_s[0] + factor * (ls_e[0] - ls_s[0])
        cpy = ls_s[1] + factor * (ls_e[1] - ls_s[1])
        cpt = math.trunc(ls_s[2] + (ls_e[2] - ls_s[2]) * factor)
        cp = [cpx, cpy, int(cpt)]
        return (cp, spherical_distance(cp, p))

    dist0 = spherical_distance(ls_s, p)
    dist1 = spherical_distance(ls_e, p)
    if dist0 < dist1:
        return (ls_s, dist0)
    return (ls_e, dist1)


def inclusion(tr1, tr2, space_th):
    """Return the sum of the distance between the two closest points of tr1 with the first and last points of tr2,
    check if tr2 is contained in tr1.
    """
    tr2_length = tr2.length()
    if tr2_length <= space_th:
        return float('infinity')

    start2 = tr2.start_point()
    end2 = tr2.end_point()

    i1_start2_point = None
    j1_end2_point = None
    i1_start2_dist = float('infinity')
    j1_end2_dist = float('infinity')

    i1 = 0
    j1 = 0

    for k in range(0, len(tr1)-1, 1):
        p1 = tr1.point_n(k)
        p2 = tr1.point_n(k+1)

        i1_start2_point_tmp = closest_point_on_segment(p1, p2, start2)
        i1_start2_dist_tmp = spherical_distance(start2, i1_start2_point_tmp)

        j1_end2_point_tmp = closest_point_on_segment(p1, p2, end2)
        j1_end2_dist_tmp = spherical_distance(end2, j1_end2_point_tmp)

        if i1_start2_dist_tmp < i1_start2_dist:
            i1_start2_dist = i1_start2_dist_tmp
            i1_start2_point = i1_start2_point_tmp
            i1 = k

        if j1_end2_dist_tmp < j1_end2_dist:
            j1_end2_dist = j1_end2_dist_tmp
            j1_end2_point = j1_end2_point_tmp
            j1 = k

    if None == i1_start2_point or None == j1_end2_point:
        return float('infinity')

    gap_i1_j1 = spherical_distance(i1_start2_point, j1_end2_point)

    if i1 >= j1 or gap_i1_j1 < space_th or (i1_start2_dist + j1_end2_dist) > tr2_length:
        return float('infinity')
    else:
        return i1_start2_dist + j1_end2_dist


def get_segments(cts, len_tr):
    segments = []
    last = cts[0] - 1
    for i in range(len(cts)):
        if cts[i] == 0:
            last = 0
            continue
        if cts[i] - last != 1:
            segments.append((last, last + 1))
            last = cts[i] - 1
        segments.append((last, cts[i]))
        last = cts[i]
    if last < len_tr-1:
        segments.append((last, last + 1))
    return segments

def has_range_intersection(cts, cte):
    cts_range = range(min(cts), max(cts))
    cte_range = range(min(cte), max(cte))
    if len(cte) == 1:
        return cte[0] in cts_range
    elif len(cts) == 1:
        return cts[0] in cte_range
    else:
        intersection = set(cts_range).intersection(set(cte_range))
        return len(intersection) != 0

def inclusion_synch(tr1, tr2, cts, cte, space_th, time_th, time_mod=86400):

    """
    cts: Closest points To Start of Tr2
    cte: Closest points To End of Tr2
    Restituisce match o None. All'interno di match, l'attributo partial_match ha il seguente significato:
            # PARTIAL MATCH 0 = IL MIGLIOR MATCH NON E' QUELLO PARZIALE / The best match is not the partial one
            # PARTIAL MATCH 1 = L'UNICO MATCH TROVATO E' QUELLO PARZIALE / The only match found is the partial one
            # PARTIAL MATCH 2 = IL MIGLIOR MATCH E' QUELLO PARZIALE / The best match is the partial one

    (*) MATCH PARZIALE = Vengono confrontati i punti iniziali e finali di tr2 con i punti iniziali e finali di tr1,
                         escludendo i punti intermedi di tr1.
    """
    #print(space_th, time_th)
    cts = sorted(cts)
    cte = sorted(cte)

    if max(cts) > min(cte) and not has_range_intersection(cts, cte):  # Two trajectories are in reverse direction: No match!
        #print("Reverse!")
        #print("Two trajectories are in reverse direction: No match!")
        return None

    start1_time = tr1.start_point()[2] % time_mod
    end1_time = tr1.end_point()[2] % time_mod

    start2_time = tr2.start_point()[2] % time_mod
    end2_time = tr2.end_point()[2] % time_mod
    #print("time ", tr1.start_point()[2], start1_time)

    # tr2 last point (plus some wasting time) is before tr1 first point
    if end2_time + time_th < start1_time:
        #print(1)
        return None

    # tr2 first point (minus some wasting time) is after tr1 last point
    if start2_time - time_th > end1_time:
        #print(2)
        return None

    # tr2 first point (plus some wasting time) is before tr1 first point and consequently
    # also before any other points of tr1
    if start2_time + time_th < start1_time:
        #print(3)
        return None

    # tr2 last point (minus some wasting time) is after tr1 last point and consequently
    # also after any other points of tr1
    if end2_time - time_th > end1_time:
        #print(4)
        return None

    tr2_length = tr2.length()
    if tr2_length <= space_th:
        #print('qui')
        #print(5)
        return None

    close_segs_start = get_segments(cts, len(tr1))
    close_segs_end = get_segments(cte, len(tr1))

    start2 = tr2.start_point()
    end2 = tr2.end_point()
    matches_with_start = []
    matches_with_end = []

    for css in close_segs_start:
        p1 = tr1.point_n(css[0])
        p2 = tr1.point_n(css[1])

        matched_point, dist_to_match = closest_point_on_segment2(p1, p2, start2)
        time_diff_to_match = int(abs((start2[2] % time_mod) - (matched_point[2] % time_mod)))

        # ToDo Must be checked later -- indexes
        matches_with_start.append((css[0], css[1], matched_point, dist_to_match, time_diff_to_match))

    for cse in close_segs_end:
        p1 = tr1.point_n(cse[0])
        p2 = tr1.point_n(cse[1])

        matched_point, dist_to_match = closest_point_on_segment2(p1, p2, end2)
        time_diff_to_match = int(abs((end2[2] % time_mod) - (matched_point[2] % time_mod)))

        # ToDo Must be checked later -- indexes
        matches_with_end.append((cse[0], cse[1], matched_point, dist_to_match, time_diff_to_match))

    #print("Matches reformed ")
    #print(matches_with_start)
    #print(matches_with_end)

    max_carpooling_index = 0
    best_match = None

    for match_s in matches_with_start:
        for match_e in matches_with_end:
            #print(match_s, match_e)
            if match_s[0] > match_e[0]:
                continue
            curr_time_diff = match_s[4] + match_e[4]
            if curr_time_diff > time_th:
                continue
            curr_sp_dist = match_s[3] + match_e[3]
            if curr_sp_dist <= space_th:
                curr_sp_dist_index, curr_time_diff_index = 1 - my_normalization(curr_sp_dist, 0, space_th), 1 - my_normalization(curr_time_diff, 0, time_th)
                curr_carpooling_index = curr_sp_dist_index * 0.5 + curr_time_diff_index * 0.5
                if max_carpooling_index < curr_carpooling_index:
                    max_carpooling_index = curr_carpooling_index
                    #print("better carpooling index")
                    #print(match_s[0]) #index
                    #print(match_e[0]) #index
                    best_match = {
                        'pickup_point': match_s[2],
                        'drop_off_point': match_e[2],
                        'dr_traj_indexes': (match_s[0], match_e[0]),
                        'space_dist_start_pickup': match_s[3],
                        'space_dist_end_drop_off': match_e[3],
                        'time_dist_start_pickup': match_s[4],
                        'time_dist_end_drop_off': match_e[4],
                        'carpooling_index': curr_carpooling_index,
                        'pass_tr_length': tr2.length()
                    }

    if best_match:
        i, j = best_match['dr_traj_indexes']
        if spherical_distance(best_match['pickup_point'], tr1.point_n(i)) > spherical_distance(best_match['pickup_point'], tr1.point_n(i+1)):
            i = i + 1
        if spherical_distance(best_match['drop_off_point'], tr1.point_n(j)) > spherical_distance(best_match['drop_off_point'], tr1.point_n(j+1)):
            j = j + 1
        indexes = (i, j)
        #print("new i, j ", indexes)
        match = best_match
        match.update({"dr_traj_indexes": indexes})
        return match
    else:
        #print("no best match")
        return None

    #print(match)



    #
    # start2 = tr2.start_point()
    # end2 = tr2.end_point()
    #
    # # con partial si intende il match solamente tra i punti iniziali e finali delle traiettorie
    # partial_dist_s, partial_dist_e, p_time_diff_s, p_time_diff_e = None, None, None, None
    # partial_dist_index, partial_time_index, partial_carpooling_index = None, None, None
    # d1 = spherical_distance(start2, tr1.start_point())
    # d2 = spherical_distance(end2, tr1.end_point())
    # td1 = abs(start1_time - start2_time)
    # td2 = abs(end1_time - end2_time)
    # if ((d1+d2) <= space_th) and ((td1+td2) <= time_th):
    #     partial_dist_s, partial_dist_e, p_time_diff_s, p_time_diff_e = d1, d2, td1, td2
    #     partial_dist_index, partial_time_index = 1 - my_normalization(d1+d2, 0, space_th), 1 - my_normalization(td1+td2, 0, time_th)
    #     partial_carpooling_index = partial_dist_index * 0.5 + partial_time_index * 0.5
    #
    # dist_s = []
    # dist_e = []
    # time_diff_s = []
    # time_diff_e = []
    # closest_point_s = []
    # closest_point_e = []
    #
    # for k in range(0, len(tr1)-1, 1):
    #     p1 = tr1.point_n(k)
    #     p2 = tr1.point_n(k+1)
    #
    #     i1_start2_point_tmp = closest_point_on_segment(p1, p2, start2)
    #     i1_start2_dist_tmp = spherical_distance(start2, i1_start2_point_tmp)
    #     i1_start2_time_diff_tmp = abs((start2[2] % time_mod) - (i1_start2_point_tmp[2] % time_mod))
    #
    #     dist_s.append(i1_start2_dist_tmp)
    #     time_diff_s.append(i1_start2_time_diff_tmp)
    #     closest_point_s.append(i1_start2_point_tmp)
    #
    #     j1_end2_point_tmp = closest_point_on_segment(p1, p2, end2)
    #     j1_end2_dist_tmp = spherical_distance(end2, j1_end2_point_tmp)
    #     j1_end2_time_diff_tmp = abs((end2[2] % time_mod) - (j1_end2_point_tmp[2] % time_mod))
    #     dist_e.append(j1_end2_dist_tmp)
    #     time_diff_e.append(j1_end2_time_diff_tmp)
    #     closest_point_e.append(j1_end2_point_tmp)
    #
    #
    # max_carpooling_index = 0
    # best_dist = float('infinity')
    # best_time_diff = None
    # indexes = (None, None)
    #
    # for i in range(len(dist_s)-1):
    #     for j in range(i+1, len(dist_s)):
    #         tmp_time_diff = time_diff_s[i] + time_diff_e[j]
    #         if tmp_time_diff > time_th:
    #             continue
    #         tmp_dist = dist_s[i] + dist_e[j]
    #         if tmp_dist <= space_th:
    #             tmp_dist_index, tmp_time_index = 1 - my_normalization(tmp_dist, 0, space_th), 1 - my_normalization(tmp_time_diff, 0, time_th)
    #             tmp_carpooling_index = tmp_dist_index * 0.5 + tmp_time_index * 0.5
    #             if tmp_carpooling_index > max_carpooling_index:
    #                 max_carpooling_index = tmp_carpooling_index
    #                 best_dist = tmp_dist
    #                 print("old code", i, j)
    #                 indexes = (i, j)
    #                 best_time_diff = tmp_time_diff
    #
    #
    # if best_dist == float('infinity'):
    #
    #     if partial_dist_s != None:
    #         # PARTIAL MATCH 1 = L'UNICO MATCH TROVATO E' QUELLO PARZIALE
    #         #print('partial match 1!')
    #         match = {
    #             'pickup_point': tr1.start_point(),
    #             'drop_off_point': tr1.end_point(),
    #             'dr_traj_indexes': (0, tr1.num_points()-1),
    #             'space_dist_start_pickup': partial_dist_s,
    #             'space_dist_end_drop_off': partial_dist_e,
    #             'time_dist_start_pickup': p_time_diff_s,
    #             'time_dist_end_drop_off': p_time_diff_e,
    #             'start_together': True,
    #             'end_together': True,
    #             'partial_match': 1,
    #             'carpooling_index': partial_carpooling_index
    #         }
    #         return match
    #     else:
    #         return None
    #
    # else:
    #     # se il match parziale è migliore tengo quello
    #     if ((partial_dist_s != None) and (max_carpooling_index < partial_carpooling_index)):
    #
    #         # PARTIAL MATCH 2 = IL MIGLIOR MATCH E' QUELLO PARZIALE
    #
    #         match = {
    #             'pickup_point': tr1.start_point(),
    #             'drop_off_point': tr1.end_point(),
    #             'dr_traj_indexes': (0, tr1.num_points()-1),
    #             'space_dist_start_pickup': partial_dist_s,
    #             'space_dist_end_drop_off': partial_dist_e,
    #             'time_dist_start_pickup': p_time_diff_s,
    #             'time_dist_end_drop_off': p_time_diff_e,
    #             'start_together': True,
    #             'end_together': True,
    #             'partial_match': 2,
    #             'carpooling_index': partial_carpooling_index
    #         }
    #         return match
    #
    #     else:
    #
    #         # PARTIAL MATCH 0 = IL MIGLIOR MATCH NON E' QUELLO PARZIALE
    #
    #         i, j = indexes
    #         pickup_point, drop_off_point = closest_point_s[i], closest_point_e[j]
    #         space_dist_start_pickup, space_dist_end_drop_off = dist_s[i], dist_e[j]
    #         time_dist_start_pickup, time_dist_end_drop_off = time_diff_s[i], time_diff_e[j]
    #
    #         newi, newj = indexes
    #         if spherical_distance(closest_point_s[i], tr1.point_n(i)) > spherical_distance(closest_point_s[i], tr1.point_n(i+1)):
    #             newi = i+1
    #         if spherical_distance(closest_point_e[j], tr1.point_n(j)) > spherical_distance(closest_point_e[j], tr1.point_n(j+1)):
    #             newj = j+1
    #         indexes = (newi, newj)
    #         print("old new indexes ", indexes)
    #
    #
    #         match = {
    #             'pickup_point': pickup_point,
    #             'drop_off_point': drop_off_point,
    #             'dr_traj_indexes': indexes,
    #             'space_dist_start_pickup': space_dist_start_pickup,
    #             'space_dist_end_drop_off': space_dist_end_drop_off,
    #             'time_dist_start_pickup': time_dist_start_pickup,
    #             'time_dist_end_drop_off': time_dist_end_drop_off,
    #             'start_together': i == 0,
    #             'end_together': j == tr1.num_points()-1,
    #             'partial_match': 0,
    #             'carpooling_index': max_carpooling_index
    #         }
    #
    #         return match
    #
    #         # return i1_start2_dist + j1_end2_dist, abs(start2_time-i1_time) + abs(end2_time-j1_time), \
    #         #        i1 == 0, j1 == tr2.num_points(), \
    #         #        i1_start2_dist, j1_end2_dist, abs(start2_time-i1_time), abs(end2_time-j1_time)







def inclusion_synch_old(tr1, tr2, space_th, time_th, time_mod=86400):

    """
    Restituisce match o None. All'interno di match, l'attributo partial_match ha il seguente significato:
            # PARTIAL MATCH 0 = IL MIGLIOR MATCH NON E' QUELLO PARZIALE
            # PARTIAL MATCH 1 = L'UNICO MATCH TROVATO E' QUELLO PARZIALE
            # PARTIAL MATCH 2 = IL MIGLIOR MATCH E' QUELLO PARZIALE

    (*) MATCH PARZIALE = Vengono confrontati i punti iniziali e finali di tr2 con i punti iniziali e finali di tr1,
                         escludendo i punti intermedi di tr1.
    """

    start1_time = tr1.start_point()[2] % time_mod
    end1_time = tr1.end_point()[2] % time_mod

    start2_time = tr2.start_point()[2] % time_mod
    end2_time = tr2.end_point()[2] % time_mod
    #print("time ", tr1.start_point()[2], start1_time)

    # tr2 last point (plus some wasting time) is before tr1 first point
    if end2_time + time_th < start1_time:
        print(1)
        return None

    # tr2 first point (minus some wasting time) is after tr1 last point
    if start2_time - time_th > end1_time:
        print(2)
        return None

    # tr2 first point (plus some wasting time) is before tr1 first point and consequently
    # also before any other points of tr1
    if start2_time + time_th < start1_time:
        print(3)
        return None

    # tr2 last point (minus some wasting time) is after tr1 last point and consequently
    # also after any other points of tr1
    if end2_time - time_th > end1_time:
        print(4)
        return None

    tr2_length = tr2.length()
    if tr2_length <= space_th:
        #print('qui')
        print(5)
        return None

    start2 = tr2.start_point()
    end2 = tr2.end_point()

    # con partial si intende il match solamente tra i punti iniziali e finali delle traiettorie
    partial_dist_s, partial_dist_e, p_time_diff_s, p_time_diff_e = None, None, None, None
    partial_dist_index, partial_time_index, partial_carpooling_index = None, None, None
    d1 = spherical_distance(start2, tr1.start_point())
    d2 = spherical_distance(end2, tr1.end_point())
    td1 = abs(start1_time - start2_time)
    td2 = abs(end1_time - end2_time)
    if ((d1+d2) <= space_th) and ((td1+td2) <= time_th):
        partial_dist_s, partial_dist_e, p_time_diff_s, p_time_diff_e = d1, d2, td1, td2
        partial_dist_index, partial_time_index = 1 - my_normalization(d1+d2, 0, space_th), 1 - my_normalization(td1+td2, 0, time_th)
        partial_carpooling_index = partial_dist_index * 0.5 + partial_time_index * 0.5

    dist_s = []
    dist_e = []
    time_diff_s = []
    time_diff_e = []
    closest_point_s = []
    closest_point_e = []

    for k in range(0, len(tr1)-1, 1):
        p1 = tr1.point_n(k)
        p2 = tr1.point_n(k+1)

        i1_start2_point_tmp = closest_point_on_segment(p1, p2, start2)
        i1_start2_dist_tmp = spherical_distance(start2, i1_start2_point_tmp)
        i1_start2_time_diff_tmp = abs((start2[2] % time_mod) - (i1_start2_point_tmp[2] % time_mod))

        dist_s.append(i1_start2_dist_tmp)
        time_diff_s.append(i1_start2_time_diff_tmp)
        closest_point_s.append(i1_start2_point_tmp)

        j1_end2_point_tmp = closest_point_on_segment(p1, p2, end2)
        j1_end2_dist_tmp = spherical_distance(end2, j1_end2_point_tmp)
        j1_end2_time_diff_tmp = abs((end2[2] % time_mod) - (j1_end2_point_tmp[2] % time_mod))
        dist_e.append(j1_end2_dist_tmp)
        time_diff_e.append(j1_end2_time_diff_tmp)
        closest_point_e.append(j1_end2_point_tmp)


    max_carpooling_index = 0
    best_dist = float('infinity')
    best_time_diff = None
    indexes = (None, None)

    for i in range(len(dist_s)-1):
        for j in range(i+1, len(dist_s)):
            tmp_time_diff = time_diff_s[i] + time_diff_e[j]
            if tmp_time_diff > time_th:
                continue
            tmp_dist = dist_s[i] + dist_e[j]
            if tmp_dist <= space_th:
                tmp_dist_index, tmp_time_index = 1 - my_normalization(tmp_dist, 0, space_th), 1 - my_normalization(tmp_time_diff, 0, time_th)
                tmp_carpooling_index = tmp_dist_index * 0.5 + tmp_time_index * 0.5
                if tmp_carpooling_index > max_carpooling_index:
                    max_carpooling_index = tmp_carpooling_index
                    best_dist = tmp_dist
                    indexes = (i, j)
                    best_time_diff = tmp_time_diff


    if best_dist == float('infinity'):

        if partial_dist_s != None:
            # PARTIAL MATCH 1 = L'UNICO MATCH TROVATO E' QUELLO PARZIALE
            #print('partial match 1!')
            match = {
                'pickup_point': tr1.start_point(),
                'drop_off_point': tr1.end_point(),
                'dr_traj_indexes': (0, tr1.num_points()-1),
                'space_dist_start_pickup': partial_dist_s,
                'space_dist_end_drop_off': partial_dist_e,
                'time_dist_start_pickup': p_time_diff_s,
                'time_dist_end_drop_off': p_time_diff_e,
                'start_together': True,
                'end_together': True,
                'partial_match': 1,
                'carpooling_index': partial_carpooling_index
            }
            return match
        else:
            return None

    else:
        # se il match parziale è migliore tengo quello
        if ((partial_dist_s != None) and (max_carpooling_index < partial_carpooling_index)):

            # PARTIAL MATCH 2 = IL MIGLIOR MATCH E' QUELLO PARZIALE

            match = {
                'pickup_point': tr1.start_point(),
                'drop_off_point': tr1.end_point(),
                'dr_traj_indexes': (0, tr1.num_points()-1),
                'space_dist_start_pickup': partial_dist_s,
                'space_dist_end_drop_off': partial_dist_e,
                'time_dist_start_pickup': p_time_diff_s,
                'time_dist_end_drop_off': p_time_diff_e,
                'start_together': True,
                'end_together': True,
                'partial_match': 2,
                'carpooling_index': partial_carpooling_index
            }
            return match

        else:

            # PARTIAL MATCH 0 = IL MIGLIOR MATCH NON E' QUELLO PARZIALE

            i, j = indexes
            pickup_point, drop_off_point = closest_point_s[i], closest_point_e[j]
            space_dist_start_pickup, space_dist_end_drop_off = dist_s[i], dist_e[j]
            time_dist_start_pickup, time_dist_end_drop_off = time_diff_s[i], time_diff_e[j]

            newi, newj = indexes
            if spherical_distance(closest_point_s[i], tr1.point_n(i)) > spherical_distance(closest_point_s[i], tr1.point_n(i+1)):
                newi = i+1
            if spherical_distance(closest_point_e[j], tr1.point_n(j)) > spherical_distance(closest_point_e[j], tr1.point_n(j+1)):
                newj = j+1
            indexes = (newi, newj)


            match = {
                'pickup_point': pickup_point,
                'drop_off_point': drop_off_point,
                'dr_traj_indexes': indexes,
                'space_dist_start_pickup': space_dist_start_pickup,
                'space_dist_end_drop_off': space_dist_end_drop_off,
                'time_dist_start_pickup': time_dist_start_pickup,
                'time_dist_end_drop_off': time_dist_end_drop_off,
                'start_together': i == 0,
                'end_together': j == tr1.num_points()-1,
                'partial_match': 0,
                'carpooling_index': max_carpooling_index
            }

            return match

            # return i1_start2_dist + j1_end2_dist, abs(start2_time-i1_time) + abs(end2_time-j1_time), \
            #        i1 == 0, j1 == tr2.num_points(), \
            #        i1_start2_dist, j1_end2_dist, abs(start2_time-i1_time), abs(end2_time-j1_time)
