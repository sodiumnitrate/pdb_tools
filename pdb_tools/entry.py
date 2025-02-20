"""
Definition of the entry class.


(entry = PDB entry in PDB lingo)
"""
import copy
import numpy as np
import gemmi
from pdb_tools.utils import get_pdbx, is_position, atom_list_to_cra_list, find_normal_to_plane_of_points
from pdb_tools.metal_site import MetalSite

import pdb

class Entry:
    """
    Definition of the Entry class that holds and manipulates data given a PDB entry.
    """
    def __init__(self, pdb_id: str):
        """
        Initializes the entry object.

        Input:
        - pdb_id (string): the 4-character PDB id
        """
        self.pdb_id = pdb_id
        
        self.structure = None
        self.fetch_structure()

        self.metal_sites = []

    def write_cif(self, file_name: str):
        """
        Function to write the .cif file corresponding to the pdb ID to file.

        Input:
        - file_name (string): file name to write structure to .cif file
        """
        self.structure.make_mmcif_document().write_file(file_name)

    def fetch_structure(self):
        """Function to fetch cif info from the PDB."""
        pdb_str = get_pdbx(self.pdb_id)
        cif_block = gemmi.cif.read_string(pdb_str)[0]
        self.structure = gemmi.make_structure_from_block(cif_block)

    def get_elements(self):
        """Function to get elements in the structure."""
        elements = []
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        el = atom.element.name
                        elements.append(el)

        return list(set(elements))

    def expand_to_P1(self):
        """Function to apply symmetry operations."""
        #TODO: how to check whether it's not NMR or any other method?
        how = gemmi.HowToNameCopiedChain.AddNumber
        self.structure.transform_to_assembly('unit_cell', how)

    def select_atoms(self, selection: str, model_idx=0):
        """
        Function to select a set of atoms, given the selection criteria.

        Input:
        - selection (str): selection string using gemmi's language
        - model_idx (int): index for model if entry has more than one model defaults to 0
        """
        if not isinstance(model_idx, int):
            raise TypeError(f"model_idx has to be an integer but {type(model_idx)} was provided.")


        # TODO: make selection language the same as VMD and convert to gemmi's language?
        sel = gemmi.Selection(selection)
        ref_atoms = []
        for i, model in enumerate(sel.models(self.structure)):
            if i != model_idx:
                continue
            for chain in sel.chains(model):
                for residue in sel.residues(chain):
                    for atom in sel.atoms(residue):
                        ref_atoms.append(atom)
        return ref_atoms

    def find_neighboring_atoms(self, ref_atom, max_dist=3, model_idx=0, tol=0.01):
        """
        Function that, given a reference atom, find neighboring atoms in a
        maximum distance of max_dist.


        Input:
        - ref_atom (gemmi.CRA): atom to search the neighborhood of
        - max_dist (float): (defaults to 3) maximum distance defining the neighborhood
        - model_idx (int): (defaults to 3) index for model if entry has more than one model 
        - tol (float): (defaults to 0.01) tolerance for distances
        """
        if not isinstance(model_idx, int):
            raise TypeError(f"model_idx has to be an integer but {type(model_idx)} was provided.")

        if isinstance(ref_atom, gemmi.CRA):
            ref_atom = ref_atom.atom
        elif not isinstance(ref_atom, gemmi.Atom):
            raise TypeError(f"ref_atom has to be a gemmi.CRA object but {type(gemmi.CRA)} was provided.")
        ns = gemmi.NeighborSearch(self.structure[model_idx],
                                  self.structure.cell,
                                  5).populate(include_h=False)

        marks = ns.find_neighbors(ref_atom, min_dist=0.01, max_dist=max_dist)

        close_atoms = []
        for mark in marks:
            cra = mark.to_cra(self.structure[model_idx])
            if cra.atom.has_altloc():
                if cra.atom.altloc == 'A':
                    close_atoms.append(cra)
            else:
                close_atoms.append(cra)

        # filter close atoms such that those that are close to each other with partial occupancy are merged
        filtered = []
        merged = []
        for i, atom1 in enumerate(close_atoms):
            if i in merged:
                continue
            if atom1.atom.occ == 1:
                filtered.append(atom1)
                continue

            occ = atom1.atom.occ
            for j, atom2 in enumerate(close_atoms):
                if occ == 1:
                    continue
                if j <= i:
                    continue
                
                if atom2.atom.occ == 1:
                    continue

                # check identity
                if atom1.atom.name != atom2.atom.name:
                    continue
                if atom1.residue.name != atom2.residue.name:
                    continue
                d_vec = atom1.atom.pos - atom2.atom.pos
                if d_vec.length() < tol:
                    occ += atom2.atom.occ
                    merged.append(j)

            atom1.atom.occ = occ
            filtered.append(atom1)

        return filtered

    def find_closest_from_each_residue(self, ref_atom, max_dist=3, model_idx=0):
        """
        Given a reference atom, find neighboring atoms such that the closest
        atom from each neighboring residue is selected.

        Input:
        - ref_atom (gemmi.CRA): reference atom
        - max_dist (float): (defaults to 3) maximum distance for neighborhood
        - model_idx (int): (defaults to 0) model index to be considered
        """
        if isinstance(ref_atom, gemmi.CRA):
            ref_atom = ref_atom.atom
        elif not isinstance(ref_atom, gemmi.Atom):
            raise TypeError
        close_atoms = self.find_neighboring_atoms(ref_atom,
                                                  max_dist=max_dist,
                                                  model_idx=model_idx)
        residues = {}
        for atom in close_atoms:
            chain = atom.chain.name
            resname = atom.residue.name
            seqid = atom.residue.seqid.num
            name = f"{chain}_{seqid}_{resname}"
            dist = self.structure.cell.find_nearest_image(ref_atom.pos, atom.atom.pos).dist()

            if name not in residues:
                residues[name] = (atom, dist)
            else:
                d_old = residues[name][1]
                if dist < d_old:
                    residues[name] = (atom, dist)

        return residues

    def find_contact_pairs(self, max_dist=4, inter_chain=False, model_idx=0):
        """
        Function to find all atom pairs that are in contact.

        Input:
        - max_dist (float): (defaults to 4)
        - inter_chain (bool): (defaults to False) whether or not contacts are to be searched on the same chain.
        - model_idx (int): (defaults to 0) model index to be considered
        """
        cs = gemmi.ContactSearch(max_dist)
        if inter_chain:
            cs.ignore = gemmi.ContactSearch.Ignore.SameChain
        else:
            cs.ignore = gemmi.ContactSearch.Ignore.AdjacentResidues

        cs.min_occupancy = 0.01

        model = self.structure[model_idx]

        # TODO: be more careful with waters
        model.remove_waters()
        ns = gemmi.NeighborSearch(model,
                                  self.structure.cell,
                                  5).populate(include_h=False)

        results = cs.find_contacts(ns)

        pairs = []
        for res in results:
            p1 = res.partner1
            p2 = res.partner2
            dist = res.dist
            pairs.append((p1, p2, dist))

        return pairs

    def find_nearest_image_distance(self, pos1, pos2):
        """
        Given two positions, find the nearest image distance between them.

        Inputs:
        - pos1: position of atom 1
        - pos2: position of atom 2
        """
        if not is_position(pos1) or not is_position(pos2):
            raise TypeError

        return self.structure.cell.find_nearest_image(gemmi.Position(pos1[0], pos1[1], pos1[2]),
                                                      gemmi.Position(pos2[0], pos2[1], pos2[2])).dist()

    def find_displacement_vector(self, pos1, pos2):
        """
        Given two positions, return the displacement vector corresponding to the
        smallest distance between the two points under PBC and symmetry operations.

        Input:
        - pos1: position of atom 1
        - pos2: position of atom 2
        """
        if not isinstance(pos1, gemmi.Position):
            pos1 = gemmi.Position(pos1[0], pos1[1], pos1[2])
        if not isinstance(pos2, gemmi.Position):
            pos1 = gemmi.Position(pos2[0], pos2[1], pos2[2])

        # TODO: is there a better way to do this??
        # TODO: you need an exception if the two atoms are the same
        i = 0
        min_dist = 0
        while True:
            try:
                dist = self.structure.cell.find_nearest_pbc_image(pos1, pos2, i)
            except IndexError:
                break
            if i == 0:
                min_dist = dist.dist()
                min_idx = i
            else:
                if dist.dist() < min_dist:
                    min_dist = dist.dist()
                    min_idx = i
            i += 1

        pos2_sym = self.structure.cell.find_nearest_pbc_position(pos1, pos2, min_idx)

        dist_vec = pos2_sym - pos1

        return [dist_vec[0], dist_vec[1], dist_vec[2]]

    def unfragmented_positions(self, atom_list):
        """
        Given a list of atoms, return a list of positions such that the nearest pbc
        image positions are returned.

        Uses the first atom in the list as a reference and brings all other atoms
        near it.

        Input:
        - atom_list (list[gemmi.CRA]): list of atoms to unfragment 
        """
        if not isinstance(atom_list, list):
            print(f"ERROR: an atom_list of {type(atom_list)} is not supported.")
            raise TypeError
        if all([isinstance(a, gemmi.CRA) for a in  atom_list]):
            positions = [a.atom.pos for a in atom_list]
        elif all([isinstance(a, gemmi.Atom) for a in atom_list]):
            positions = [a.pos for a in atom_list]
        elif all([isinstance(a, gemmi.Atom) or isinstance(a, gemmi.CRA) for a in atom_list]):
            positions = []
            for a in atom_list:
                if isinstance(a, gemmi.Atom):
                    positions.append(a.pos)
                elif isinstance(a, gemmi.CRA):
                    positions.append(a.atom.pos)
        else:
            raise TypeError

        ref_pos = positions[0]
        new_positions = [np.array([ref_pos[0], ref_pos[1], ref_pos[2]])]
        for pos in positions[1:]:
            dist_vec = self.find_displacement_vector(ref_pos, pos)
            new_pos = np.array(dist_vec) + np.array([ref_pos[0], ref_pos[1], ref_pos[2]])
            new_positions.append(np.copy(new_pos))

        return np.array(new_positions)

    def find_center_of_mass(self, atom_list=None, selection_string=None, unweighted=False, model_idx=0):
        """
        Find center of mass of a given atom list or selection of atoms.

        If atom_list is provided, selection_string and model_idx are rendered irrelevant.

        Input:
        - atom_list (list[gemmi.CRA]): list containing atoms
        - selection_string (str): string that defines atom selection using gemmi's language
        - unweighted (bool): (defaults to False) if True, treats weights of all atoms to be the same
        - model_idx (int): (defaults to 0) model index for structure to be considered
        """
        if atom_list is None and selection_string is None:
            print("ERROR: atom_list and selection can't be both None.")
            raise ValueError

        if atom_list is None:
            atom_list = self.select_atoms(selection_string, model_idx=model_idx)

        positions = self.unfragmented_positions(atom_list)
        com = np.array([0., 0., 0.])
        tot_w = 0
        for i, atom in enumerate(atom_list):
            if isinstance(atom, gemmi.CRA):
                w = atom.atom.element.weight
            elif isinstance(atom, gemmi.Atom):
                w = atom.element.weight
            else:
                raise TypeError
            if unweighted:
                w = 1
            tot_w += w
            com += positions[i] * w

        com /= tot_w
        return com

    def get_metal_sites(self, metal_element: str, model_idx=0):
        """
        Given the metal element, find and compile a list of metal sites.

        Skips H and C.

        Inputs:
        - metal_element (str): symbol for element to be considered
        - model_idx (int): (defaults to 0) model index for structure to be considered
        """
        elements = self.get_elements()
        if metal_element not in elements:
            print(f"ERROR: {metal_element} not in the structure with PDB ID: {self.pdb_id} and model {model_idx}.")
            raise ValueError

        self.metal_sites = []

        selection_string = f"[{metal_element}]"
        atoms = self.select_atoms(selection_string, model_idx=model_idx)

        for atom in atoms:
            atom_cra = atom_list_to_cra_list([atom], self.structure, model_idx=model_idx)
            neigs = self.find_neighboring_atoms(atom, max_dist=3, model_idx=model_idx)
            atom_list = [atom]
            for neig in neigs:
                if neig.atom.pos == atom_cra[0].atom.pos and neig.atom.name == atom_cra[0].atom.name:
                    continue
                if neig.atom.element.name != 'H' and neig.atom.element.name != 'C':
                    atom_list.append(neig)

            positions = self.unfragmented_positions(atom_list)

            new_metal = MetalSite(atom,
                                  atom_list[1:],
                                  positions)
            self.metal_sites.append(new_metal)
    
    def get_nearby_aromatic(self, ref_atom, max_dist=7, model_idx=0, tol=0.01):
        """
        Given a reference atom, find nearby aromatic rings
        that belong to Y, F, W, or H.

        Input:
        - ref_atom (gemmi.CRA): reference atom
        - max_dist (float): (defaults to 7) maximum distance to consider for neighborhood
        - model_idx (int): (defaults to 0): model index for structure to be considered
        - tol (float): (defaults to 0.01) distance tolerance
        """
        atom_list = self.find_neighboring_atoms(ref_atom, max_dist=max_dist, model_idx=model_idx, tol=tol)
        aromatic_residues = []
        for atom in atom_list:
            if atom.residue.name in ['TRP', 'TYR', 'HIS', 'PHE']:
                residue = (atom.residue.name,
                           atom.chain.name,
                           atom.residue.seqid.num)
                if residue not in aromatic_residues:
                    aromatic_residues.append(residue)

        return aromatic_residues

    def get_aromatic_ring_center_and_normal(self, aromatic_residue, model_idx=0, trp_ring=0):
        """
        Given an aromatic residue, find its aromatic ring center(s) and
        the normal to the plane containing the ring(s).

        Input:
        - aromatic_residue = (str: resname, str: chain, int: seqid)

        Output:
        - com: center of mass vector
        - normal: normal vector

        TODO: add functionality to deal with a gemmi.Residue object.
        """
        resname = aromatic_residue[0]
        if resname not in ['TYR','PHE','HIS','TRP']:
            print("ERROR: the residue selected isn't aromatic.")
            raise ValueError

        chain = aromatic_residue[1]
        seqid = aromatic_residue[2]

        rings = {'TYR':['CG', 'CD2', 'CE2', 'CZ', 'CE1', 'CD1'],
                 'PHE':['CG', 'CD2', 'CE2', 'CZ', 'CE1', 'CD1'],
                 'HIS':['CG', 'ND1', 'CE1', 'NE2', 'CD2'],
                 'TRP_0':['CD1', 'CG', 'CD2', 'CE2', 'NE1'],
                 'TRP_1':['CD2', 'CE3', 'CZ3', 'CH2', 'CZ2', 'CE2']}

        if resname == 'TRP':
            resname += f"_{trp_ring}"

        selection_string = f"/{model_idx+1}/{chain}/{seqid}/"

        selected_atoms = []
        for name in rings[resname]:
            s = selection_string + f"{name}[{name[0]}]"
            atom = self.select_atoms(s)
            selected_atoms += atom

        positions = self.unfragmented_positions(selected_atoms)
        normal = find_normal_to_plane_of_points(positions)
        com = self.find_center_of_mass(selected_atoms, unweighted=True)

        return com, normal

    def find_aromatics_near_metal_sites(self, cutoff=7,
                                              max_d0=5.5, max_theta0=55,
                                              model_idx=0, tol=0.01):
        """
        After finding metal sites, find aromatic residues around them.

        Input:
        - cutoff (float): maximum distance to consider
        - max_d0 (float): maximum distance to consider between metal and aromatic ring
        - max_theta (float): maximum angle (degrees) to consider between aromatic ring plane and metal
        - model_idx (int): (defaults to 0) 
        - tol (float): (defaults to 0.01)
        """
        if len(self.metal_sites) == 0:
            print("ERROR: couldn't find any metal sites. Did you run get_metal_sites() first?")
            raise ValueError

        for site in self.metal_sites:
            metal_atom = site.metal
            aromatic_residues = self.get_nearby_aromatic(metal_atom, max_dist=cutoff,
                                                         model_idx=model_idx,
                                                         tol=tol)
            nearby = []
            d0_vals = []
            theta0_vals = []
            for res in aromatic_residues:
                if res[0] == "TRP":
                    n = 2
                else:
                    n = 1
                for i in range(n):   
                    com, normal = self.get_aromatic_ring_center_and_normal(res, model_idx=model_idx, trp_ring=i)
                    disp_vec = self.find_displacement_vector(com, metal_atom.pos)
                    d0 = np.linalg.norm(disp_vec)
                    theta0_up = np.arccos(np.dot(disp_vec, normal) / (d0 * np.linalg.norm(normal))) * 180 / np.pi
                    theta0_down = np.arccos(np.dot(disp_vec, -1*normal) / (d0 * np.linalg.norm(normal))) * 180 / np.pi
                    theta0 = min(theta0_up, theta0_down)
                    if d0 <= max_d0 and theta0 <= max_theta0:
                        nearby.append(res)
                        d0_vals.append(d0)
                        theta0_vals.append(theta0)

            site.nearby_aromatic = copy.copy(nearby)
            site.d0 = copy.copy(d0_vals)
            site.theta0 = copy.copy(theta0_vals)
