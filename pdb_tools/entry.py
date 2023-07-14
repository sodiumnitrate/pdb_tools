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
        how = gemmi.HowToNameCopiedChain.AddNumber
        self.structure.transform_to_assembly('unit_cell', how)
