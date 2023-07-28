"""
Simple implementation of ICP, based on https://github.com/ClayFlannigan/icp

(Simple is enough, as the max number of points we have is 10. ICP is even an
overkill, but I wanted to learn about it.)
"""
import numpy as np
from itertools import permutations
from sklearn.neighbors import NearestNeighbors

from pdb_tools.utils import is_two_d_array_of_floats

import pdb

def nearest_neighbor(src, dst):
    """
    Given two sets of points, find the nearest neighbor distance
    and index of nearest neighbor for each point.
    """

    neigh = NearestNeighbors(n_neighbors=1)
    neigh.fit(dst)
    distances, indices = neigh.kneighbors(src, return_distance=True)
    return distances.ravel(), indices.ravel()

def nearest_neighbor_best_match(pts1, pts2):
    """
    Given two sets of points, find matching that minimizes distance.
    """
    if pts1.shape != pts2.shape:
        print("ERROR: input point sets must have the same dimensions.")
        raise ValueError

    dist_dict = {}
    for i in range(len(pts1)):
        pos_i = pts1[i,:]
        for j in range(i, len(pts1)):
            pos_j = pts2[j,:]
            d = np.linalg.norm(pos_i - pos_j)
            dist_dict[(i,j)] = d

    min_dist = np.inf
    min_matching = None
    for matching in permutations(list(range(len(pts1))), len(pts1)):
        indices = list(matching)
        distances = []
        for i, idx in enumerate(indices):
            key = (min(i,idx), max(i,idx))
            distances.append(dist_dict[key])

        distances = np.array(distances)
        rmsd = np.mean(distances**2)
        if rmsd < min_dist:
            min_dist = rmsd
            min_matching = indices

    return min_dist, min_matching


def best_fit_no_translation(points_1, points_2):
    """
    Find the best fit rotation that maps points_1 to points2.
    """
    n_dim = points_1.shape[1]

    # rotation matrix
    H = np.dot(points_1.T, points_2)
    U, S, Vt = np.linalg.svd(H)
    R = np.dot(Vt.T, U.T)

    # special reflection case
    if np.linalg.det(R) < 0:
        Vt[n_dim-1,:] *= -1
        R = np.dot(Vt.T, U.T)

    return R

def icp_no_translation(points_1, points_2, max_iter=20, tol=0.0001):
    """
    Given a set of points_1, and a set of points_2, find the best-fit
    transform mapping points_1 to points_2.

    Input Nxn_dim array.
    """

    # check data types
    if not is_two_d_array_of_floats(points_1) or not is_two_d_array_of_floats(points_2):
        raise TypeError
    points_1 = np.array(points_1)
    points_2 = np.array(points_2)

    if points_1.shape != points_2.shape:
        raise ValueError

    n_dim = points_1.shape[1]
    if n_dim != 3:
        print("ERROR: data isn't 3d.")
        raise ValueError

    src = np.copy(points_1.T)
    dst = np.copy(points_2.T)

    prev_error = 0

    for i in range(max_iter):
        # find the nearest neighbors
        rmsd, indices = nearest_neighbor_best_match(src.T, dst.T)

        # best fit transformation
        R = best_fit_no_translation(src.T, dst[:,indices].T)

        # update current source
        src = np.dot(R, src)

        print(rmsd, indices)

        # check error
        if np.abs(prev_error - rmsd) < tol:
            break

        prev_error = rmsd
    
    # final transformation
    R = best_fit_no_translation(points_1, src.T)

    return R, rmsd, i

    
def best_fit_transform(points_1, points_2):
    """
    Find the best fit transform that maps points_1 to points_2.
    """
    n_dim = points_1.shape[1]

    # translate to centroids
    com_1 = np.mean(points_1, axis=0)
    com_2 = np.mean(points_2, axis=0)
    points_1p = points_1 - com_1
    points_2p = points_2 - com_2

    # rotation matrix
    H = np.dot(points_1p.T, points_2p)
    U, S, Vt = np.linalg.svd(H)
    R = np.dot(Vt.T, U.T)

    # special reflection case
    if np.linalg.det(R) < 0:
        Vt[n_dim-1,:] *= -1
        R = np.dot(Vt.T, U.T)

    # translation
    t = com_2.T - np.dot(R, com_1.T)

    # homogeneous transformation
    T = np.identity(n_dim + 1)
    T[:n_dim, :n_dim] = R
    T[:n_dim, n_dim] = t

    return T, R, t


def iterative_closest_point(points_1, points_2, max_iter=20, tol=0.0001):
    """
    Given a set of points_1, and a set of points_2, find the best-fit
    transform mapping points_1 to points_2.

    Input Nxn_dim array.
    """

    # check data types
    if not is_two_d_array_of_floats(points_1) or not is_two_d_array_of_floats(points_2):
        raise TypeError
    points_1 = np.array(points_1)
    points_2 = np.array(points_2)

    if points_1.shape != points_2.shape:
        raise ValueError

    n_dim = points_1.shape[1]
    if n_dim != 3:
        print("ERROR: data isn't 3d.")
        raise ValueError

    # expand coord arrays so that you can apply the entire transformation in one go
    src = np.ones((n_dim + 1, points_1.shape[0]))
    dst = np.ones((n_dim + 1, points_2.shape[0]))
    src[:n_dim, :] = np.copy(points_1.T)
    dst[:n_dim, :] = np.copy(points_2.T)

    prev_error = 0

    for i in range(max_iter):
        # find the nearest neighbors
        distances, indices = nearest_neighbor(src[:n_dim,:].T, dst[:n_dim,:].T)

        # best fit transformation
        T, _, _ = best_fit_transform(src[:n_dim,:].T, dst[:n_dim,indices].T)

        # update current source
        src = np.dot(T, src)

        # check error
        mean_error = np.sqrt(np.mean(distances**2))
        if np.abs(prev_error - mean_error) < tol:
            break

        prev_error = mean_error
    
    # final transformation
    T, R, t = best_fit_transform(points_1, src[:n_dim, :].T)

    return R, distances, i, t
