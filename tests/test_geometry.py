"""
Unit tests for the functions and entities in geometry.py.
"""
from pdb_tools.geometry import *
import numpy as np

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
        min_geom, min_rmsd, distorted = find_coordination_geometry(ligands, center)
        assert min_geom == 'tet'
        assert min_rmsd < 1e-10
        assert not distorted

    def test_find_coordination_geometry_2(self):
        ligands = trigonal_prism()
        center = [0, 0, 0]
        min_geom, min_rmsd, distorted = find_coordination_geometry(ligands, center)

        assert min_geom == 'tpr'
        assert min_rmsd < 1e-10

    def test_compare_two_geometries_cn2(self):
        rmsd = compare_two_geometries('lin', 'lin')
        assert rmsd == 0
        rmsd = compare_two_geometries('lin', 'trv')
        assert np.isclose(rmsd*3, 1.553, atol=0.001)

        rmsd2 = compare_two_geometries('trv', 'lin')
        assert np.isclose(rmsd, rmsd2)

    def test_compare_two_geometries_cn3(self):
        cn = 3
        rmsds = {}
        for g1 in coordination_numbers[cn]:
            for g2 in coordination_numbers[cn]:
                rmsds[(g1,g2)]= compare_two_geometries(g1, g2)

        for key, val in rmsds.items():
            print(key, val)
            if key[0] == key[1]:
                assert np.isclose(val, 0)

            rmsd2 = rmsds[(key[1], key[0])]
            assert np.isclose(rmsd2, val)

    def test_compare_two_geometries_cn4(self):
        cn = 4
        rmsds = {}
        for g1 in coordination_numbers[cn]:
            for g2 in coordination_numbers[cn]:
                rmsds[(g1,g2)] = compare_two_geometries(g1, g2)

        for key, val in rmsds.items():
            print(key, val)
            if key[0] == key[1]:
                assert np.isclose(val, 0)

            rmsd2 = rmsds[(key[1], key[0])]
            assert np.isclose(rmsd2, val)

    def test_compare_two_geometries_cn5(self):
        cn = 5
        rmsds = {}
        for g1 in coordination_numbers[cn]:
            for g2 in coordination_numbers[cn]:
                rmsds[(g1,g2)] = compare_two_geometries(g1, g2)

        for key, val in rmsds.items():
            print(key, val)
            if key[0] == key[1]:
                assert np.isclose(val, 0)

            rmsd2 = rmsds[(key[1], key[0])]
            assert np.isclose(rmsd2, val)

    def test_1(self):
        g1 = 'spl'
        g2 = 'bva'

        pts1 = np.array(ligand_coords[g1], dtype=np.float_)
        pts2 = np.array(ligand_coords[g2], dtype=np.float_)

        rmsd, best_rot, pts2p = best_match(pts1, pts2)

        rmsd2, best_rot2, pts1p = best_match(pts2, pts1)
        assert np.isclose(rmsd, rmsd2)

        rmsd, best_rot, pts2p = best_match(pts1, pts1)
        assert rmsd == 0