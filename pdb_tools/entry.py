"""
Definition of the entry class.

(entry = PDB entry in PDB lingo)
"""
import gemmi
from pdb_tools.utils import get_pdbx

class Entry:
    def __init__(self, pdb_id):
        self.pdb_id = pdb_id
        
        self.structure = None
        self.fetch_structure()

    def write_cif(self, file_name):
        """Function to write the .cif file corresponding to the .pdb ID to file."""
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
        close_atoms = self.find_neighboring_atoms(ref_atom,
                                                  max_dist=max_dist,
                                                  model_idx=model_idx)
        residues = {}
        for atom in close_atoms:
            chain = atom.chain.name
            resname = atom.residue.name
            seqid = atom.residue.seqid.num
            name = f"{chain}_{seqid}_{resname}"
            dist = self.structure.cell.find_nearest_pbc_image(ref_atom.pos, atom.atom.pos, 0).dist()

            if name not in residues:
                residues[name] = (atom, dist)
            else:
                d_old = residues[name][1]
                if dist < d_old:
                    residues[name] = (atom, dist)

        return residues
