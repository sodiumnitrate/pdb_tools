"""
Unit tests for the ICP implementation in icp.py.
"""
import numpy as np

import pdb

from pdb_tools.icp import iterative_closest_point, best_fit_transform
from pdb_tools.geometry import tetrahedral
from scipy.spatial.transform import Rotation as R

class TestICP:
    def test_best_fit_transform(test):
        points_1 = np.array(tetrahedral()).T
        rot = R.from_euler('z', 90, degrees=True).as_matrix()
        points_2 = rot.dot(points_1)

        T, Rmatrix, t = best_fit_transform(points_1.T, points_2.T)
        assert len(Rmatrix) == 3
        assert len(t) == 3

        assert np.abs(np.sum(rot - Rmatrix)) < 1e-10

        assert np.abs(np.sum(t)) < 1e-10


    def test_icp_1(self):
        points_1 = np.array(tetrahedral()).T
        rot = R.from_euler('z', 90, degrees=True).as_matrix()
        points_2 = rot.dot(points_1)

        points_1 = points_1.T
        points_2 = points_2.T

        rotation, distances, n_iter, t = iterative_closest_point(points_1, points_2)

        assert np.abs(np.sum(t)) < 1e-10
        assert np.abs(np.sum(distances)) < 1e-10

        assert np.abs(np.abs(np.linalg.det(rotation)) - 1 ) < 1e-10