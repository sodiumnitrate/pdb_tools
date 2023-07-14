"""
Unit tests for utils.
"""

from pdb_tools.utils import *
from pdb_tools.entry import Entry

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