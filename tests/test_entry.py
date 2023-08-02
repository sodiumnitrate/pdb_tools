"""
Unit tests for the Entry class.
"""
from math import sqrt
import numpy as np
from pdb_tools.entry import Entry
from pdb_tools.utils import cra_from_atom
import gemmi

import pdb

class TestEntry:
    def test_init(self):
        entry = Entry('6q6b')

        assert entry.structure 

    def test_get_elements(self):
        entry = Entry('6q6b')
        elements = entry.get_elements()
        assert len(elements) > 0
        assert 'C' in elements

    def test_write_cif(self, tmpdir):
        entry = Entry('6q6b')
        file = tmpdir.mkdir("sub").join('test.cif')
        entry.write_cif(file.strpath)
        assert len(file.read()) > 0


    def test_expand_to_P1(self):
        entry = Entry('6q6b')
        n_old = entry.structure[0].count_atom_sites()
        entry.expand_to_P1()
        n_new = entry.structure[0].count_atom_sites()

        assert n_new > n_old

    def test_select_atoms(self):
        entry = Entry('6q6b')
        atoms = entry.select_atoms('(CU1)')

        assert len(atoms) == 59
        assert atoms[0].name == 'CU'

    def test_find_neighboring_atoms(self):
        entry = Entry('6q6b')
        ref_atom = entry.structure[0][0][0][0]
        neighboring_atoms = entry.find_neighboring_atoms(ref_atom)

        assert len(neighboring_atoms) == 5

    def test_find_neighboring_atoms_2(self):
        entry = Entry('6q6b')
        ref_atom = entry.structure[0][0][0][0]
        ref_atom = cra_from_atom(ref_atom, entry.structure)
        neighboring_atoms = entry.find_neighboring_atoms(ref_atom)

        assert len(neighboring_atoms) == 5

    def test_closest_from_each_residue(self):
        entry = Entry('6q6b')
        ref_atom = entry.structure[0][0][0][0]
        
        residues = entry.find_closest_from_each_residue(ref_atom)
        assert len(residues) == 2

    def test_closest_from_each_residue(self):
        entry = Entry('6q6b')
        ref_atom = entry.structure[0][0][0][0]
        ref_atom = cra_from_atom(ref_atom, entry.structure)
        residues = entry.find_closest_from_each_residue(ref_atom)
        assert len(residues) == 2

    def test_find_contact_pairs(self):
        entry = Entry('6q6b')
        results = entry.find_contact_pairs()

        assert len(results) == 7246
        assert isinstance(results[0], tuple)


    def test_closest_2(self):
        pdb_id = '2c9p'
        entry = Entry(pdb_id)
        entry2 = Entry(pdb_id)
        entry2.expand_to_P1()
        
        coppers = entry.select_atoms('(CU)')
        coppers2 = entry2.select_atoms('(CU)')
        
        for i in range(len(coppers)):
            neigs = entry.find_closest_from_each_residue(coppers[i])
            neigs2 = entry2.find_closest_from_each_residue(coppers2[i])
            
            assert len(neigs) == len(neigs2)

    def test_find_displacement_vector(self):
        entry = Entry('2c9p')

        coppers = entry.select_atoms('(CU)')
        copper = coppers[1]

        pos1 = copper.pos
        neigs = entry.find_closest_from_each_residue(copper)
        pos2 = neigs['A_48_HIS'][0].atom.pos

        dist = entry.structure.cell.find_nearest_image(pos1, pos2).dist()

        dist_vec = entry.find_displacement_vector(pos1, pos2)

        assert abs(dist - sqrt(dist_vec[0]**2 + dist_vec[1] **2 + dist_vec[2] ** 2)) < 1e-10

    def test_unfragmented_positions(self):
        entry = Entry('2c9p')

        coppers = entry.select_atoms('(CU)')
        copper = coppers[1]

        neigs = entry.find_neighboring_atoms(copper)
        positions = entry.unfragmented_positions([copper] + neigs)

        for i in range(len(positions)):
            pos_i = positions[i]
            for j in range(i+1, len(positions)):
                pos_j = positions[j]
                dist_1 = np.linalg.norm(pos_i - pos_j)

                dist_2 = entry.find_nearest_image_distance(pos_i, pos_j)

                assert abs(dist_1 - dist_2) < 1e-10

    def test_find_com(self):
        entry = Entry('2c9p')

        coppers = entry.select_atoms('(CU)')
        copper = coppers[1]

        neigs = entry.find_neighboring_atoms(copper)

        com = entry.find_center_of_mass(neigs)

        assert np.linalg.norm(np.array([copper.pos[0], copper.pos[1], copper.pos[2]]) - com) < 1

    def test_get_metal_sites(self):
        entry = Entry('2c9p')
        entry.get_metal_sites('Cu')

        assert len(entry.metal_sites) > 0

        for site in entry.metal_sites:
            for l in site.ligands:
                assert l.atom.element.name != 'C'
                assert l.atom.element.name != 'H'

    def test_get_metal_sites_2(self):
        entry = Entry('4DYX')
        atoms = entry.select_atoms('[Cu]')
        cu_1 = atoms[0]
        cu_2 = atoms[1]
        neigs_2 = entry.find_neighboring_atoms(cu_2)

        assert len(neigs_2) == 6

    def test_get_nearby_aromatic(self):
        entry = Entry('1GST')
        atoms = entry.select_atoms('/1/A/55/OD1[O]')
        aromatic_residues = entry.get_nearby_aromatic(atoms[0])

        assert ('PHE', 'B', 140) in aromatic_residues
        assert ('PHE', 'A', 56) in aromatic_residues

    def test_get_aromatic_ring(self):
        entry = Entry('1GST')
        res = ('PHE', 'A', 56)
        com, normal = entry.get_aromatic_ring_center_and_normal(res)

        assert len(com) == 3
        assert len(normal) == 3

    def test_find_aromatic_near_metal(self):
        entry = Entry('1LLA')
        entry.get_metal_sites('Cu')
        entry.find_aromatics_near_metal_sites(max_d0=5.5, max_theta0=55)

        for site in entry.metal_sites:
            assert site.nearby_aromatic is not None
            assert site.d0 is not None
            assert site.theta0 is not None

        assert len(entry.metal_sites[0].nearby_aromatic) == 0
        assert len(entry.metal_sites[1].nearby_aromatic) == 2

    def test_find_aromatic_near_metal_2(self):
        entry = Entry('2gli')
        entry.get_metal_sites('Co')
        entry.find_aromatics_near_metal_sites(max_d0=5.5, max_theta0=55)

        for site in entry.metal_sites:
            assert site.nearby_aromatic is not None
            assert site.d0 is not None
            assert site.theta0 is not None

        assert len(entry.metal_sites[0].nearby_aromatic) == 1
        assert entry.metal_sites[0].nearby_aromatic[0] == ('TRP', 'A', 108)