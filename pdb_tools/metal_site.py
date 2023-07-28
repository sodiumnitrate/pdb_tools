"""
Functions and tools to find metal site given a pdb structure.
"""

from pdb_tools.utils import get_pdbx
from pdb_tools.geometry import *

class MetalSite:
    def __init__(self, metal, ligands, unfragmented_coords):
        self.metal = metal
        self.ligands = ligands
        self.unfragmented_coords = unfragmented_coords
        self.coordination_number = len(ligands)

        self.assigned_geometry = None
        self.rmsd = None

    def assign_geometry(self):
        ligands = self.unfragmented_coords[1:]
        center = self.unfragmented_coords[0]

        min_geom, min_rmsd, distorted = find_coordination_geometry(ligands, center)

        self.assigned_geometry = min_geom
        self.rmsd = min_rmsd
        self.distorted = distorted

    def print_geometry(self):
        if self.assigned_geometry is None:
            self.assign_geometry()
        
        if self.assign_geometry is None:
            print("No geometry has been assigned.")
            return

        name = coordination_geometry_names[self.assign_geometry]
        if self.distorted:
            name += " (distorted)"
        name += f" rmsd={self.rmsd}"

        print(name)

    def plot_metal_site(self):
        plot_geometry(self.unfragmented_coords, plot_zero=False)

    def check_all_geometries(self):
        ligands = self.unfragmented_coords[1:]
        center = self.unfragmented_coords[0]

        result = []
        coordination_number = len(ligands)
        for geom in coordination_numbers[coordination_number]:
            rmsd, l_AT, l_DT = check_coordination_geometry(ligands, center, geom)
            result.append((geom, rmsd, l_AT, l_DT))

        return result