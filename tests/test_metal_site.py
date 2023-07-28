"""
Unit tests for the MetalSite class.
"""

from pdb_tools.entry import Entry

class TestMetalSite:
    def test_assign_geometry(self):
        entry = Entry('2c9p')
        entry.get_metal_sites('Cu')
        for site in entry.metal_sites:
            site.assign_geometry()