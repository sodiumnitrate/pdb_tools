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

        # whether there's an aromatic residue nearby
        self.nearby_aromatic = None
        self.aromatic_d0 = None
        self.aromatic_theta0 = None

    def assign_geometry(self):
        ligands = self.unfragmented_coords[1:]
        center = self.unfragmented_coords[0]

        min_geom, min_rmsd, distorted = find_coordination_geometry(ligands, center)

        self.assigned_geometry = min_geom
        self.rmsd = min_rmsd
        self.distorted = distorted

        if self.assigned_geometry is None:
            self.assigned_geometry = 'irr'

    def print_geometry(self):
        if self.assigned_geometry is None:
            self.assign_geometry()

        name = coordination_geometry_names[self.assigned_geometry]
        if self.distorted is not None:
            if self.distorted:
                name += " (distorted)"
            else:
                name += " (regular)"
        if self.rmsd is not None:
            name += f" rmsd={self.rmsd:.3f}"

        return name

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