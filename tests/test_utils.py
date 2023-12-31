"""
Unit tests for utils.
"""
from pdb_tools.utils import *
from pdb_tools.entry import Entry

import pdb

class TestUtils:
    def test_get_pdbx(self):
        ans = get_pdbx("6q6b")
        assert isinstance(ans, bytes)
        assert len(ans) > 0

    def test_chimera_selection_from_atom_list(self):
        entry = Entry('6q6b')
        ref_atom = entry.structure[0][0][0][0]

        neighboring_atoms = entry.find_neighboring_atoms(ref_atom)

        res_string = chimera_selection_from_atom_list(neighboring_atoms)
        assert isinstance(res_string, str)

        #TODO: how to test validity?

    def test_vmd_selection_from_atom_list(self):
        entry = Entry('6q6b')
        ref_atom = entry.structure[0][0][0][0]

        neighboring_atoms = entry.find_neighboring_atoms(ref_atom)

        res_string = vmd_selection_from_atom_list(neighboring_atoms)
        assert isinstance(res_string, str)

        #TODO: how to test validity?

    def test_cra_from_atom(self):
        entry = Entry('6q6b')
        ref_atom = entry.structure[0][0][0][0]
        assert isinstance(ref_atom, gemmi.Atom)

        cra = cra_from_atom(ref_atom, entry.structure)

        assert cra.atom == ref_atom

    def test_atom_list_to_cra_list(self):
        entry = Entry('6q6b')
        atom_list = [entry.structure[0][0][0][1], entry.structure[0][0][0][5]]

        cra_list = atom_list_to_cra_list(atom_list, structure=entry.structure)

        assert cra_list[0].atom == atom_list[0]

    def test_vmd_selection_from_atom_list_2(self):
        entry = Entry('6q6b')
        atom_list = [entry.structure[0][0][0][1], entry.structure[0][0][0][5]]

        res_string = vmd_selection_from_atom_list(atom_list, structure=entry.structure)

        assert isinstance(res_string, str)

    def test_are_points_roughly_planar(self):
        points = [[0,0,0], [-5,0,0.1], [0,5,0], [5,0,-0.1], [-5, -5, 0.5]]

        assert are_points_roughly_planar(points)

    def test_is_position(self):
        a = [0, 1, 20]
        assert is_position(a)

        a = np.array([0, 1, 10])
        assert is_position(a)

        a = gemmi.Position(0, 10, 20)
        assert is_position(a)

    def test_is_two_d_array_of_floats(self):
        array = [[0.1,0.5,1.0], [0.5, 1.0, 1.5]]
        assert is_two_d_array_of_floats(array)

        array = np.array(array)
        assert is_two_d_array_of_floats(array)

        array = [[[5.1,2.2],[0.1,0.2]],[[0.1,2.0],[11., 5.]]]
        assert not is_two_d_array_of_floats(array)

        array = np.array(array)
        assert not is_two_d_array_of_floats(array)

        array = np.array([2,3.])
        assert not is_two_d_array_of_floats(array)

        array = [2.1, 0.3]
        assert not is_two_d_array_of_floats(array)

        array = [['foo','fofo'],['bar','baba']]
        assert not is_two_d_array_of_floats(array)

        array = [[2,0.2],[0,0]]
        assert not is_two_d_array_of_floats(array)

    def test_find_normal(self):
        points = [[0,0,0], [-5,0,0.1], [0,5,0], [5,0,-0.1], [-5, -5, 0.5]]
        normal = find_normal_to_plane_of_points(points)
        assert len(normal) == 3
        assert isinstance(normal, np.ndarray)