"""
Definition of the entry class.

(entry = PDB entry in PDB lingo)
"""
import numpy as np
import gemmi
from pdb_tools.utils import get_pdbx, is_position

class Entry:
    def __init__(self, pdb_id):
        self.pdb_id = pdb_id
        
        self.structure = None
        self.fetch_structure()

    def write_cif(self, file_name):
        """Function to write the .cif file corresponding to the pdb ID to file."""
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

    def select_atoms(self, selection, model_idx=0):
        """Function to select a set of atoms, given the selection criteria."""
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

    def find_neighboring_atoms(self, ref_atom, max_dist=3, model_idx=0):
        """
        Function that, given a reference atom, find neighboring atoms in a
        maximum distance of max_dist.
        """
        if isinstance(ref_atom, gemmi.CRA):
            ref_atom = ref_atom.atom
        elif not isinstance(ref_atom, gemmi.Atom):
            raise TypeError
        ns = gemmi.NeighborSearch(self.structure[model_idx],
                                  self.structure.cell,
                                  5).populate(include_h=False)

        marks = ns.find_neighbors(ref_atom, min_dist=0.01, max_dist=max_dist)

        close_atoms = []
        for mark in marks:
            cra = mark.to_cra(self.structure[model_idx])
            close_atoms.append(cra)

        return close_atoms

    def find_closest_from_each_residue(self, ref_atom, max_dist=3, model_idx=0):
        """
        Given a reference atom, find neighboring atoms such that the closest
        atom from each neighboring residue is selected.
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
        """Function to find all atom pairs that are in contact."""
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
        """
        if not is_position(pos1) or not is_position(pos2):
            raise TypeError

        return self.structure.cell.find_nearest_image(gemmi.Position(pos1[0], pos1[1], pos1[2]),
                                                      gemmi.Position(pos2[0], pos2[1], pos2[2])).dist()

    def find_displacement_vector(self, pos1, pos2):
        """
        Given two positions, return the displacement vector corresponding to the
        smallest distance between the two points under PBC and symmetry operations.
        """

        # TODO: is there a better way to do this??
        i = 0
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
        Given a list of atom, return a list of positions such that the nearest pbc
        image positions are returned.

        Uses the first atom in the list as a reference and brings all other atoms
        near it.
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

    def find_center_of_mass(self, atom_list=None, selection_string=None, model_idx=0):
        """
        Find center of mass of a given atom list or selection of atoms.

        If atom_list is provided, selection_string and model_idx are rendered irrelevant.
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
            tot_w += w
            com += positions[i] * w

        com /= tot_w
        return com