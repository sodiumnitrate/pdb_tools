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