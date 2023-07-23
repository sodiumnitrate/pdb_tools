"""
Unit tests for the functions and entities in geometry.py.
"""
from pdb_tools.geometry import *
import numpy as np

import pdb

class TestGeometry:
    def test_labels(self):
        labels_from_CN = []
        for val in coordination_numbers.values():
            labels_from_CN += val

        assert coordination_geometry_labels == labels_from_CN

        labels_from_coord_dict = list(ligand_coords.keys())
        assert set(labels_from_coord_dict) == set(coordination_geometry_labels)

    def test_coordination_numbers(self):
        for lbl, coords in ligand_coords.items():
            cn = label_to_CN[lbl]
            assert len(coords) == cn

    def test_unit_vectors(self):
        for lbl, coords in ligand_coords.items():
            for vec in coords:
                l = np.linalg.norm(vec)
                assert np.isclose(l, 1)

    def test_min_angle(self):
        for lbl, coords in ligand_coords.items():
            print(lbl)
            for i in range(len(coords)):
                for j in range(i+1, len(coords)):
                    v1 = np.array(coords[i])
                    v2 = np.array(coords[j])
                    d = v1.dot(v2)
                    if np.isclose(d, -1):
                        d = -1
                    elif np.isclose(d, 1):
                        d = 1
                    theta = np.arccos(d)
                    assert theta > np.pi/18


    def test_find_coordination_geometry_1(self):
        ligands = tetrahedral()
        center = [0, 0, 0]

        min_geom, min_rmsd = find_coordination_geometry(ligands, center)
        assert min_geom == 'tet'
        assert min_rmsd < 1e-10

    def test_find_coordination_geometry_2(self):
        ligands = tetrahedral()
        center = [0.1, 0.1, 0.1]

        min_geom, min_rmsd = find_coordination_geometry(ligands, center)
        assert min_geom == 'tet'
        assert min_rmsd > 1e-10

    def test_find_coordination_geometry_3(self):
        ligands = np.array(tetrahedral()) * 3
        center = [0.1, 0.1, 0.1]

        min_geom, min_rmsd = find_coordination_geometry(ligands, center)
        assert min_geom == 'tet'

    def test_find_coordination_geometry_4(self):
        ligands = trigonal_prism()
        center = [0, 0, 0]
        min_geom, min_rmsd = find_coordination_geometry(ligands, center)

        assert min_geom == 'tpr'
        assert min_rmsd < 1e-10