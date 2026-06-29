import os
import sys
import shutil
import json
import numpy as np
import pandas as pd
import networkx as nx
import osmnx as ox
import geopandas as gpd

from tqdm import tqdm
import multiprocessing.pool as mpp
from sklearn.preprocessing import MaxAbsScaler

sys.path.append('..')
sys.path.append('../libs')
from personalized_routing import (_worker_compute_user_gaps,
                                  _init_global_variables,
                                  _istarmap,
                                  _log,
                                  attributes)

from trajectory_analysis import _worker_collect_user_attributes


def learn_user_pref_allones(G: nx.MultiDiGraph,
                            routes: list,
                            cost_attributes: list,
                            alphas_path: str,
                            real_path: str,
                            results_dir: str,
                            iterations: int,
                            num_processes: int = None,
                            verbose: bool = False):
    """
    Same as the function in personalized_routing.py, but does not distinguish
    between clusters of users, and initializes all betas to 1.0 for all cost
    attributes.
    """
    # TODO replace paths with actual dataframes/dicts

    # Creating folder for results
    if os.path.exists(results_dir):
        _log(f"Warning: results directory {results_dir} already exists. It will be overwritten.", verbose)
        shutil.rmtree(results_dir)
    os.makedirs(results_dir)
    
    # Conversion to indices/matrix format for scipy
    nodes = list(G.nodes())
    n_nodes = len(G.nodes)
    node_to_idx = {node: i for i, node in enumerate(nodes)}
    idx_to_node = {i: node for i, node in enumerate(nodes)}

    edges = list(G.edges(data=True))
    u_indices = np.array([node_to_idx[u] for u, _, _ in edges])
    v_indices = np.array([node_to_idx[v] for _, v, _ in edges])

    attribute_values = []
    for attr in attributes:
        attribute_values.append(np.array([d[attr] for _, _, d in edges]))
    attribute_values = np.array(attribute_values)

    base_cost = np.nanmedian(attribute_values[attributes.index('travel_time')])
    scaler = MaxAbsScaler()
    attribute_values = scaler.fit_transform(attribute_values)

    real_means = pd.read_csv(real_path)[['uid'] + [f'{attr}_mean' for attr in attributes]]
    real_means.columns = ['uid'] + attributes

    alphas = json.load(open(alphas_path, 'r'))

    # Patching multiprocessing to allow starmap + tqdm progress bar
    mpp.Pool.istarmap = _istarmap
    pool = mpp.Pool(processes=num_processes,
                    initializer=_init_global_variables,
                    initargs=(G,
                              n_nodes,
                              node_to_idx,
                              idx_to_node,
                              u_indices,
                              v_indices,
                              attribute_values,
                              cost_attributes,
                              base_cost)
    )

    betas = {
        cost: 1.0 for cost in cost_attributes
    }

    # List of MAE values across all sub-iterations
    mae = []
    best_maes = {attr: np.inf for attr in cost_attributes}

    for i in range(iterations * len(cost_attributes)):
        c_idx = i % len(cost_attributes) # index in cost_attributes
        a_idx = attributes.index(cost_attributes[c_idx]) # index in attributes
        
        _log(f"It {i // len(cost_attributes)}, adjusting attribute {cost_attributes[c_idx]}", verbose)

        c_beta = betas[cost_attributes[c_idx]]
        beta_gaps = {} # Collects gaps for all users for each candidate beta
        beta_overlaps = {} # Collects overlaps for all users for each candidate beta

        for candidate_beta in [c_beta/10,
                            c_beta/2,
                            c_beta,
                            c_beta*10/2,
                            c_beta*10]:
            _log(f"Testing beta {candidate_beta}", verbose)

            user_results = list(
                tqdm(pool.istarmap(_worker_compute_user_gaps,
                                [(uid,
                                    user_routes,
                                    real_means[real_means['uid'] == uid][attributes].values[0],
                                    candidate_beta,
                                    c_idx,
                                    i,
                                    alphas[str(uid)],
                                    betas) 
                                    for uid, user_routes in routes.items()]),
                    desc='Users',
                    total=len(routes),
                    disable=not verbose)
            )
            beta_gaps[candidate_beta] = []
            beta_overlaps[candidate_beta] = []
            for gap, overlaps in user_results:
                beta_gaps[candidate_beta].append(gap)
                beta_overlaps[candidate_beta].extend(overlaps)

        # Checking if any beta improves the MAE for the current cost attribute
        best_gaps = None
        best_overlap = None
        for b, g in beta_gaps.items():
            g = np.array(g)
            b_mae = np.nanmean(np.abs(g[:,a_idx]))

            if b_mae < best_maes[cost_attributes[c_idx]]:
                # Improvement
                best_maes[cost_attributes[c_idx]] = b_mae
                betas[cost_attributes[c_idx]] = b
                best_gaps = g
                best_overlap = beta_overlaps[b]

        if best_gaps is None:
            # No improvement found
            _log(f"No improvement found for attribute  {cost_attributes[c_idx]}", verbose)
        else:
            _log(f"Found improvement for attribute {cost_attributes[c_idx]}: {betas[cost_attributes[c_idx]]}, MAE: {best_maes[cost_attributes[c_idx]]}", verbose)

            # Saving intermediate relative gaps
            pd.DataFrame(best_gaps, columns=attributes).to_csv(os.path.join(results_dir, f'gaps_iter_{i}.csv'), index=False)
        
            # Saving intermediate overlaps
            pd.DataFrame(best_overlap, columns=['uid', 'tid', 'overlap']).to_csv(os.path.join(results_dir, f'overlaps_iter_{i}.csv'), index=False)
            _log('Mean Overlap: {:.4f}'.format(np.mean(best_overlap, axis=0)[2]), verbose)
            
            # Saving current beta values
            json.dump(betas, open(os.path.join(results_dir, f'betas_iter_{i}.json'), 'w'))

        mae.append(best_maes[cost_attributes[c_idx]])

    _log(f"Learned betas: {betas}", verbose)

    json.dump(betas, open(os.path.join(results_dir, f'final_betas.json'), 'w'))
    np.save(os.path.join(results_dir, f'mae.npy'), np.array(mae))

    pool.close()
    pool.join()


def test_betas(G: nx.MultiDiGraph,
               routes: dict,
               cost_attributes: list,
               ratios: dict,
               betas: dict,
               real_path: str,
               results_path: str,
               num_processes: int = None,
               verbose: bool = False):
    """
    Test a set of beta values for the cost function and compute the relative
    gaps between the real route attributes and synthetic route attributes. The
    relative gaps are written to file.
    """

    # Conversion to matrix format 
    nodes = list(G.nodes())
    n_nodes = len(G.nodes)
    node_to_idx = {node: i for i, node in enumerate(nodes)}
    idx_to_node = {i: node for i, node in enumerate(nodes)}

    edges = list(G.edges(data=True))
    u_indices = np.array([node_to_idx[u] for u, _, _ in edges])
    v_indices = np.array([node_to_idx[v] for _, v, _ in edges])

    attribute_values = []
    for attr in attributes:
        attribute_values.append(np.array([d[attr] for _, _, d in edges]))
    attribute_values = np.array(attribute_values)

    base_cost = np.nanmedian(attribute_values[attributes.index('travel_time')])
    scaler = MaxAbsScaler()
    attribute_values = scaler.fit_transform(attribute_values)

    mpp.Pool.istarmap = _istarmap
    pool = mpp.Pool(processes=num_processes,
                    initializer=_init_global_variables,
                    initargs=(G,
                              n_nodes,
                              node_to_idx,
                              idx_to_node,
                              u_indices,
                              v_indices,
                              attribute_values,
                              cost_attributes,
                              base_cost)
    )

    real_means = pd.read_csv(real_path)[['uid'] + [f'{attr}_mean' for attr in attributes]]
    real_means.columns = ['uid'] + attributes

    with open(results_path, 'w') as f:
        f.write(f'uid,{",".join(attributes)}\n')
        for uid, result in zip(routes.keys(), list(tqdm(pool.istarmap(_worker_compute_user_gaps,
                                                                      [(uid,
                                                                        user_routes,
                                                                        real_means[real_means['uid'] == uid][attributes].values[0],
                                                                        None,
                                                                        None,
                                                                        4,
                                                                        ratios[str(uid)],
                                                                        betas) for uid, user_routes in routes.items()]),
                                                        desc='Users',
                                                        total=len(routes),
                                                        disable=not verbose))):
            gap, _ = result
            if gap is not None:
                f.write(f'{uid},' + ','.join([str(g) for g in gap]) + '\n')

    pool.close()
    pool.join()


def random_search_betas(G: nx.MultiDiGraph,
                        routes: dict,
                        path_results: str,
                        filename_json: str,
                        num_processes: int = None,
                        verbose: bool = False):
    """
    Perform a random search over a space of beta values for the cost function.
    """

    beta_values = [0.1, 0.5, 0.25, 1.0, 2.5, 5.0, 10.0, 25.0, 50.0, 100.0]
    cost_attributes = ['residential', 'landuse_industrial', 'landuse_commercial', 'flow_low', 'water']

    # Generate 200 random combinations
    combinations = set()
    while len(combinations) < 200:
        beta_combination = tuple(np.random.choice(beta_values, len(cost_attributes)))
        if os.path.exists(os.path.join(path_results, f'gap_beta_{"_".join(map(str, beta_combination))}.csv')):
            continue
        combinations.add(beta_combination)

    # Convert graph to matrix format
    nodes = list(G.nodes())
    node_to_idx = {node: i for i, node in enumerate(nodes)}
    idx_to_node = {i: node for i, node in enumerate(nodes)}

    edges = list(G.edges(data=True))
    u_indices = np.array([node_to_idx[u] for u, _, _ in edges])
    v_indices = np.array([node_to_idx[v] for _, v, _ in edges])

    attribute_values = []
    for attr in attributes:
        attribute_values.append(np.array([d[attr] for _, _, d in edges]))
    attribute_values = np.array(attribute_values)

    base_cost = np.nanmedian(attribute_values[attributes.index('travel_time')])

    scaler = MaxAbsScaler()
    attribute_values = scaler.fit_transform(attribute_values)

    mpp.Pool.istarmap = _istarmap
    pool = mpp.Pool(processes=num_processes,
                    initializer=_init_global_variables,
                    initargs=(G,
                              len(G.nodes),
                              node_to_idx,
                              idx_to_node,
                              u_indices,
                              v_indices,
                              attribute_values,
                              cost_attributes,
                              base_cost)
    )

    ratios = json.load(open(filename_json, 'r'))
    
    if not os.path.exists(path_results):
        os.makedirs(path_results)

    df_real = pd.read_csv('../statistics/edge_attributes_real.csv.gz', compression='gzip')
    real_means = df_real.groupby(['uid'], sort=False)[attributes].mean().reset_index()

    for beta_combination in tqdm.tqdm(combinations, desc='Beta combinations', disable=not verbose):
        user_attributes = list(tqdm.tqdm(
            pool.istarmap(
                _worker_collect_user_attributes,
                [(uid,
                user_routes,
                ratios[str(uid)],
                beta_combination) for uid, user_routes in routes.items()]),
            total=len(routes),
            desc='Users')
        )

        gaps = []
        for uid, user_values in user_attributes:
            if uid not in real_means['uid'].values:
                continue

            xbar_real = real_means[real_means['uid'] == uid][attributes].values[0]
            xbar_combo = np.array([np.nanmean(user_values[attr]) for attr in attributes])

            gap_u = np.where(
                xbar_combo != 0.0,
                (xbar_real - xbar_combo) / xbar_combo,
                np.where(
                    xbar_real != 0.0,
                    np.nan,
                    0.0
                )
            )
            gaps.append(np.append(uid, gap_u))
        gaps = np.array(gaps)
        pd.DataFrame(gaps, columns=['uid'] + attributes).to_csv(os.path.join(path_results, f'gap_beta_{"_".join(map(str, beta_combination))}.csv'), index=False)

    pool.close()
    pool.join()
