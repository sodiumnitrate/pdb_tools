"""
This file contains a bunch of utilities for the module.
"""
from pypdb.clients.pdb import pdb_client
import gemmi

def get_pdbx(pdb_id):
    return pdb_client.get_pdb_file(pdb_id, pdb_client.PDBFileType.CIF, True)

def chimera_selection_from_atom_list(atom_list):
    """
    Given a list of gemmi CRA objects, generate selection
    string that works with chimera.
    """
    #TODO: atom to CRA conversion on the fly?

    selection = []
    for atom in atom_list:
        chain = atom.chain.name
        resname = atom.residue.name
        seqid = atom.residue.seqid.num
        a = atom.atom.name

        atom_selection = f"(/{chain} & :{seqid} & @{a})"
        selection.append(atom_selection)

    return "select " + " | ".join(selection)

def vmd_selection_from_atom_list(atom_list):
    """
    Given a list of gemmi CRA objects, generate selection
    string that works with VMD.
    """
    #TODO: atom to CRA conversion on the fly?

    selection = []
    for atom in atom_list:
        chain = atom.chain.name
        resname = atom.residue.name
        seqid = atom.residue.seqid.num
        a = atom.atom.name

        atom_selection = f"(resid {seqid} and name {a} and chain {chain})"
        selection.append(atom_selection)

    return " or ".join(selection)