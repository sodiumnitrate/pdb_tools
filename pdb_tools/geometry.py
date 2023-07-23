"""
Functions and tools to do the necessary computational geometry.
"""
import copy
import numpy as np
import matplotlib.pyplot as plt

from pdb_tools.icp import iterative_closest_point

# coordination labels from Andreini et al.
coordination_geometry_labels = ['lin', 'trv', 'tri', 'tev', 'spv', 'tet', 'spl',
                                'bva', 'bvp', 'pyv', 'spy', 'tbp', 'tpv', 'oct',
                                'tpr', 'pva', 'pvp', 'cof', 'con', 'ctf', 'ctn',
                                'pbp', 'coc', 'ctp', 'hva', 'hvp', 'cuv', 'sav',
                                'hbp', 'cub', 'sqa', 'boc', 'bts', 'btt', 'ttp',
                                'csa']

CN_2 = ['lin', 'trv']
CN_3 = ['tri', 'tev', 'spv']
CN_4 = ['tet', 'spl', 'bva', 'bvp', 'pyv']
CN_5 = ['spy', 'tbp', 'tpv']
CN_6 = ['oct', 'tpr', 'pva', 'pvp', 'cof', 'con', 'ctf', 'ctn']
CN_7 = ['pbp', 'coc', 'ctp', 'hva', 'hvp', 'cuv', 'sav']
CN_8 = ['hbp', 'cub', 'sqa', 'boc', 'bts', 'btt']
CN_9 = ['ttp', 'csa']

coordination_numbers = {2: CN_2, 3: CN_3, 4: CN_4, 5: CN_5, 6: CN_6, 7: CN_7, 8: CN_8, 9: CN_9}
label_to_CN = {}
for key, value in coordination_numbers.items():
    for el in value:
        label_to_CN[el] = key

ligand_coords = {'lin': [[-1., 0, 0], [1, 0, 0]],
                 'trv': [[-1*np.sqrt(3)/2, -0.5, 0], [1, 0, 0]],
                 'tri': [[-1*np.sqrt(3)/2, -0.5, 0], [1, 0, 0], [-1*np.sqrt(3)/2, 0.5, 0]]}

def tetrahedral():
    # function to calculate and return unit vectors for tetrahedral geom
    vec1 = [1, 0, 0]
    cos_theta = -1/3
    theta = np.arccos(cos_theta)
    vec2 = [np.cos(theta), 0, np.sin(theta)]
    v3z = (np.cos(theta)-np.cos(theta)**2)/np.sin(theta)
    v3y = np.sqrt(1-np.cos(theta)**2-v3z**2)
    vec3 = [np.cos(theta), v3y, v3z]
    vec4 = [np.cos(theta), -v3y, v3z]
    return [vec1, vec2, vec3, vec4]

tet = tetrahedral()
ligand_coords['tev'] = copy.copy(tet[:-1])
ligand_coords['tet'] = copy.copy(tet)

# square planars
ligand_coords['spv'] = [[1.,0.,0.], [-1., 0., 0.], [0., -1., 0.]]
ligand_coords['spl'] = [[1.,0.,0.], [-1., 0., 0.], [0., -1., 0.], [0., 1., 0.]]

# trigonal bipyramidals
trigonal_bipyramidal = [[0,0,1.], [0,0,-1], [1,0,0],
                        [-1*np.sqrt(3)/2, -0.5, 0],
                        [-1*np.sqrt(3)/2, 0.5, 0]]

ligand_coords['tbp'] = copy.copy(trigonal_bipyramidal)
ligand_coords['bva'] = copy.copy(trigonal_bipyramidal[1:])
ligand_coords['bvp'] = copy.copy(trigonal_bipyramidal[:-1])

# square pyramidals
sq_py = [[1.,0.,0.], [-1.,0.,0.], [0.,1.,0.],[0.,-1.,0.],[0.,0.,1.]]
ligand_coords['spy'] = copy.copy(sq_py)
ligand_coords['pyv'] = copy.copy(sq_py[1:])

# tpv/tpr (trigonal prism)
def trigonal_prism():
    l = np.sqrt(12/13)
    h = l / np.sqrt(3)

    v1 = [0, l, 0.5*h]
    v2 = [l*np.sqrt(3)/2, -l/2, 0.5*h]
    v3 = [-l*np.sqrt(3)/2, -l/2, 0.5*h]
    v4 = [0, l, -0.5*h]
    v5 = [l*np.sqrt(3)/2, -l/2, -0.5*h]
    v6 = [-l*np.sqrt(3)/2, -l/2, -0.5*h]

    return [v1, v2, v3, v4, v5, v6]

tpr = trigonal_prism()
ligand_coords['tpr'] = copy.copy(tpr)
ligand_coords['tpv'] = copy.copy(tpr[:-1])

# octahedral
ligand_coords['oct'] = copy.copy(sq_py) + [[0,0,-1]]

# pbp/pva/pvp (pentagonal bipyramid)
angle = 2 * np.pi / 5
pbp = [[1,0,0],
       [0,0,1],
       [0,0,-1],
       [np.cos(angle), np.sin(angle), 0],
       [-np.cos(angle/2), np.sin(angle/2), 0],
       [-np.cos(angle/2), -np.sin(angle/2), 0],
       [np.cos(angle), -np.sin(angle), 0]]
ligand_coords['pbp'] = copy.copy(pbp)
ligand_coords['pva'] = copy.copy(pbp[:1] + pbp[2:])
ligand_coords['pvp'] = copy.copy(pbp[:-1])

# coc/cof/con (octahedron, face mono-capped)
ligand_coords['coc'] = copy.copy(ligand_coords['oct']) + [[1/np.sqrt(3), 1/np.sqrt(3), 1/np.sqrt(3)]]
ligand_coords['cof'] = copy.copy(ligand_coords['coc'][1:])
ligand_coords['con'] = copy.copy(ligand_coords['oct'][:-1]) + [[1/np.sqrt(3), 1/np.sqrt(3), 1/np.sqrt(3)]]

# hbp/hva/hvp (hexagonal bipyramid)
ligand_coords['hbp'] = [[0,0,1], [0,0,-1],
                        [1,0,0], [-1,0,0],
                        [0.5, np.sqrt(3)/2, 0],
                        [-0.5, np.sqrt(3)/2, 0],
                        [-0.5, -np.sqrt(3)/2, 0],
                        [0.5, -np.sqrt(3)/2, 0]]
ligand_coords['hva'] = copy.copy(ligand_coords['hbp'][1:])
ligand_coords['hvp'] = copy.copy(ligand_coords['hbp'][:-1])

# cub/cuv (cube)
l = 1/np.sqrt(3)
cub = [[l, l, l], [-l, l, l], [-l, -l, l], [l, -l, l],
       [l, l, -l], [-l, l, -l], [-l, -l, -l], [l, -l, -l]]
ligand_coords['cub'] = copy.copy(cub)
ligand_coords['cuv'] = copy.copy(cub[:-1])

# sqa/sav (square antiprism) --> twist the top of cub 45°
l = 1/np.sqrt(3)
ll = 2/np.sqrt(6)
sqa = [[ll,0,l],
       [-ll,0,l],
       [0,ll,l],
       [0,-ll,l],
       [l,l,-l], [-l,l,-l], [-l,-l,-l], [l,-l,-l]]
ligand_coords['sqa'] = copy.copy(sqa)
ligand_coords['sav'] = copy.copy(sqa[1:])

# boc (octahedron, trans-bicapped)
ligand_coords['boc'] = copy.copy(ligand_coords['coc']) + [[-1/np.sqrt(3), -1/np.sqrt(3), -1/np.sqrt(3)]]

# btt (trigonal prism, triangular-face bicapped) --> add an extra atom to all the triangular faces
ligand_coords['btt'] = copy.copy(ligand_coords['tpr']) + [[0,0,1], [0,0,-1]]

# ttp (trigonal prism, square-face tricapped) --> add an extra atom to all the square faces
ligand_coords['ttp'] = copy.copy(ligand_coords['tpr']) + [[0, -1, 0], [np.sqrt(3)/2, 0.5, 0], [-np.sqrt(3)/2, 0.5,0]]

# bts (trigonal prism, square-face bicapped) --> add an extra atom to two of the square faces
ligand_coords['bts'] = copy.copy(ligand_coords['tpr']) + [[0, -1, 0], [np.sqrt(3)/2, 0.5, 0]]

# ctp/ctf/ctn (trigonal prism, square-face monocapped)
ligand_coords['ctp'] = copy.copy(ligand_coords['tpr']) + [[0, -1, 0]]
ligand_coords['ctn'] = copy.copy(ligand_coords['ctp'][1:])
ligand_coords['ctf'] = copy.copy(ligand_coords['ctp'][:1]) + copy.copy(ligand_coords['ctp'][2:])

# csa (square antiprism, square-face monocapped) --> add an extra atom to the top
ligand_coords['csa'] = copy.copy(sqa) + [[0,0,1]]

def find_coordination_geometry(ligands, center):
    """
    Given an array of coordinates for the ligands (N by 3), and given
    the coordinates of the center (1 by 3), determine the geometry that
    most closely matches it.
    """
    # get coordination number
    coordination_number = len(ligands)

    coords = normalize_bond_lengths(ligands, center)

    possibilities = []
    for geom in coordination_numbers[coordination_number]:
        geom_coords = np.array(ligand_coords[geom] + [[0,0,0]], dtype=np.float_)

        _, distances, _ = iterative_closest_point(coords, geom_coords)
        rmsd = np.sqrt(np.sum(distances**2)) / (coordination_number + 1)

        possibilities.append((geom, rmsd))

    min_rmsd = possibilities[0][1]
    min_geom = possibilities[0][0]
    for geom, rmsd in possibilities:
        if rmsd < min_rmsd:
            min_rmsd = rmsd
            min_geom = geom

    return min_geom, min_rmsd

def normalize_bond_lengths(ligands, center):
    """
    Given an array of coordinates for the ligands (N by 3), and given
    the coordinates of the center (1 by 3), move center to [0,0,0],
    and normalize bond lengths.
    """
    coords = []
    for l in ligands:
        v1 = np.array(l) - np.array(center, dtype=np.float_)
        v1 /= np.linalg.norm(v1)
        coords.append(v1)

    coords.append(np.array([0.,0.,0.]))

    return np.array(coords)

def plot_geometry(coords):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.set_aspect('equal')
    for coord in coords:
        ax.scatter3D(coord[0], coord[1], coord[2], color='blue')

    ax.scatter3D([0], [0], [0], color='red')
    plt.show()