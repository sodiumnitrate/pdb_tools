"""
This file contains a bunch of utilities for the module.
"""
from pypdb.clients.pdb import pdb_client

def get_pdbx(pdb_id):
    return pdb_client.get_pdb_file(pdb_id, pdb_client.PDBFileType.CIF, True)
