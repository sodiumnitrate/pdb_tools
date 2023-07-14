"""
Unit tests for the Entry class.
"""


from pdb_tools.entry import Entry


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

