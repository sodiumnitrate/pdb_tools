"""
Simple implementation of ICP, based on https://github.com/ClayFlannigan/icp

(Simple is enough, as the max number of points we have is 10. ICP is even an
overkill, but I wanted to learn about it.)
"""
import numpy as np
from sklearn.neighbors import NearestNeighbors

from pdb_tools.utils import is_two_d_array_of_floats

def nearest_neighbor(src, dst):
    """
    Given two sets of points, find the nearest neighbor distance
    and index of nearest neighbor for each point.
    """

    neigh = NearestNeighbors(n_neighbors=1)
    neigh.fit(dst)
    distances, indices = neigh.kneighbors(src, return_distance=True)
    return distances.ravel(), indices.ravel()

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


def iterative_closest_point(points_1, points_2, max_iter=20, tol=0.001):
    """
    Given a set of points_1, and a set of points_2, find the best-fit
    transform mapping points_1 to points_2.

    Given Nxn_dim or n_dimxN arrays, the assumption is that N>n_dim.
    """

    # check data types
    if not is_two_d_array_of_floats(points_1) or not is_two_d_array_of_floats(points_2):
        raise TypeError
    points_1 = np.array(points_1)
    points_2 = np.array(points_2)

    shape_1 = points_1.shape
    shape_2 = points_2.shape

    if shape_1[0] < shape_1[1]:
        points_1 = points_1.T

    if shape_2[0] < shape_2[1]:
        points_2 = points_2.T

    if points_1.shape != points_2.shape:
        raise ValueError

    n_dim = points_1.shape[1]

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
        mean_error = np.mean(distances)
        if np.abs(prev_error - mean_error) < tol:
            break

        prev_error = mean_error
    
    # final transformation
    T, _, _ = best_fit_transform(points_1, src[:n_dim, :].T)

    return T, distances, i
