from pdb_tools.utils import *

class TestUtils:
    def test_get_pdbx(self):
        ans = get_pdbx("6q6b")
        assert isinstance(ans, bytes)
        assert len(ans) > 0
