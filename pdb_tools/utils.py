"""
This file contains a bunch of utilities for the module.
"""
from pypdb.clients.pdb import pdb_client
import numpy as np
import gemmi
from skspatial.objects import Plane, Points

def get_pdbx(pdb_id):
    return pdb_client.get_pdb_file(pdb_id, pdb_client.PDBFileType.CIF, True)

def chimera_selection_from_atom_list(atom_list, structure=None, model_idx=0, tol=0.0001):
    """
    Given a list of gemmi CRA objects, generate selection
    string that works with chimera.
    """
    atom_list = atom_list_to_cra_list(atom_list, structure=structure, model_idx=model_idx, tol=tol)

    selection = []
    for atom in atom_list:
        chain = atom.chain.name
        resname = atom.residue.name
        seqid = atom.residue.seqid.num
        a = atom.atom.name

        atom_selection = f"(/{chain} & :{seqid} & @{a})"
        selection.append(atom_selection)

    return "select " + " | ".join(selection)

def atom_list_to_cra_list(atom_list, structure=None, model_idx=0, tol=0.0001):
    """
    Function to check atom list.
    """
    if isinstance(atom_list, list):
        if isinstance(atom_list[0], gemmi.Atom):
            if structure is None:
                print("ERROR: a list of gemmi.Atom objects provided, but no structure object is provided.")
                raise ValueError
        elif not isinstance(atom_list[0], gemmi.CRA):
            print(f"ERROR: wrong input type {type(atom_list[0])}.")
            raise TypeError
    elif isinstance(atom_list, gemmi.Atom) or isinstance(atom_list, gemmi.CRA):
        atom_list = [atom_list]
    else:
        print(f"ERROR: unknown type {type(atom_list)}.")
        raise TypeError
    
    return [cra_from_atom(atom, structure, model_idx=model_idx, tol=tol) for atom in atom_list]

def vmd_selection_from_atom_list(atom_list, structure=None, model_idx=0, tol=0.0001):
    """
    Given a list of gemmi CRA objects, generate selection
    string that works with VMD.
    """
    atom_list = atom_list_to_cra_list(atom_list, structure=structure, model_idx=model_idx, tol=tol)
    
    selection = []
    for atom in atom_list:
        chain = atom.chain.name
        seqid = atom.residue.seqid.num
        a = atom.atom.name

        atom_selection = f"(resid {seqid} and name {a} and chain {chain})"
        selection.append(atom_selection)

    return " or ".join(selection)

def cra_from_atom(atom, structure, model_idx=0, tol=0.0001):
    """
    Given a gemmi model, and a gemmi atom object, find gemmi CRA object.
    """
    if isinstance(atom, gemmi.CRA):
        return atom

    pos = atom.pos

    return cra_from_atom_position(pos,
                                  structure,
                                  model_idx=model_idx,
                                  tol=tol)

def cra_from_atom_position(pos, structure, model_idx=0, tol=0.0001):
    """
    Given an atom position, find gemmi CRA object.
    """
    if isinstance(pos, list):
        if len(pos) != 3:
            raise ValueError
        pos = gemmi.Position(pos[0], pos[1], pos[2])
    elif isinstance(pos, np.ndarray):
        if len(pos) != 3:
            raise ValueError
        pos = gemmi.Position(pos[0], pos[1], pos[2])
    elif not isinstance(pos, gemmi.Position):
        raise ValueError

    ns = gemmi.NeighborSearch(structure[model_idx],
                              structure.cell,
                              5).populate()

    marks = ns.find_atoms(pos, radius=tol)

    if len(marks) > 1:
        print("WARNING: more than one atom found. Try reducing the tolerance.")
    elif len(marks) == 0:
        print("WARNING: no atoms found. Try increasing the tolerance.")
        return None

    return marks[0].to_cra(structure[model_idx])

def are_points_roughly_planar(points, max_err=0.1):
    """
    Given a set of points, fit a plane to them, and compute distances
    to the plane to figure out if they are roughly planar.

    https://stackoverflow.com/a/72384583
    """

    #TODO: check type and size of points
    points = Points(points)
    best_fit_plane = Plane.best_fit(points)

    distances = [best_fit_plane.distance_point(point) for point in points]
    error = np.sqrt(np.dot(distances, distances) / len(distances))

    if error > max_err:
        return False

    return True

def is_position(vec):
    """
    Given a vector, find out if it's consistent with a position vector
    in 3d.
    """

    if isinstance(vec, gemmi.Position):
        return True

    if isinstance(vec, list) or isinstance(vec, np.ndarray):
        if all([isinstance(a, (int, float, np.int_, np.float_)) for a in vec]):
            if len(vec) == 3:
                return True

    return False