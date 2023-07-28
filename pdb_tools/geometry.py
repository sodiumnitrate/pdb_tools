"""
Functions and tools to do the necessary computational geometry.

TODO:
- find a way to optimize geometry finding. Currently, CN>7 is a bit of a problem.
"""

import copy
import numpy as np
import matplotlib.pyplot as plt
from itertools import permutations

#from pdb_tools.icp import icp_no_translation

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

coordination_geometry_names = {label:'' for label in coordination_geometry_labels}
coordination_geometry_names['lin'] = 'Linear'
coordination_geometry_names['trv'] = 'Trigonal plane with valence'
coordination_geometry_names['tri'] = 'Trigonal plane'
coordination_geometry_names['tev'] = 'Tetrahedron with a vacancy'
coordination_geometry_names['spv'] = 'Square plane with a vacancy'
coordination_geometry_names['tet'] = 'Tetrahedron'
coordination_geometry_names['spl'] = 'Square plane'
coordination_geometry_names['bva'] = 'Trigonal bipyramid with a vacancy (axial)'
coordination_geometry_names['bvp'] = 'Trigonal bipyramid with a vacancy (equatorial)'
coordination_geometry_names['pyv'] = 'Square pyramid with a vacancy (equatorial)'
coordination_geometry_names['spy'] = 'Square pyramid'
coordination_geometry_names['tbp'] = 'Trigonal bipyramid'
coordination_geometry_names['tpv'] = 'Trigonal prism with a vacancy'
coordination_geometry_names['oct'] = 'Octahedron'
coordination_geometry_names['tpr'] = 'Trigonal prism'
coordination_geometry_names['pva'] = 'Pentagonal bipyramid with a vacancy (axial)'
coordination_geometry_names['pvp'] = 'Pentagonal bipyramid with a vacancy (equatorial)'
coordination_geometry_names['cof'] = 'Octahedron, face monocapped with a vacancy (capped face)'
coordination_geometry_names['con'] = 'Octahedron, face monocapped with a vacancy (non-capped face)'
coordination_geometry_names['ctf'] = 'Trigonal prism, square-face monocapped with a vacancy (capped face)'
coordination_geometry_names['ctn'] = 'Trigonal prism, square-face monocapped with a vacancy (non-capped edge)'
coordination_geometry_names['pbp'] = 'Pentagonal bipyramid'
coordination_geometry_names['coc'] = 'Octahedron, face monocapped'
coordination_geometry_names['ctp'] = 'Trigonal prism, square-face monocapped'
coordination_geometry_names['hva'] = 'Hexagonal bipyramid with a vacancy (axial)'
coordination_geometry_names['hvp'] = 'Hexagonal bipyramid with a vacancy (equatorial)'
coordination_geometry_names['cuv'] = 'Cube with a vacancy'
coordination_geometry_names['sav'] = 'Square antiprism with a vacancy'
coordination_geometry_names['hbp'] = 'Hexagonal bipyramid'
coordination_geometry_names['cub'] = 'Cube'
coordination_geometry_names['sqa'] = 'Square antiprism'
coordination_geometry_names['boc'] = 'Octahedron, trans-bicapped'
coordination_geometry_names['bts'] = 'Trigonal prism, square-face bicapped'
coordination_geometry_names['btt'] = 'Trigonal prism, triangular-face bicapped'
coordination_geometry_names['ttp'] = 'Trigonal prism, square-face tricapped'
coordination_geometry_names['csa'] = 'Square antiprism, square-face monocapped'

coordination_numbers = {2: CN_2, 3: CN_3, 4: CN_4, 5: CN_5, 6: CN_6, 7: CN_7, 8: CN_8, 9: CN_9}
label_to_CN = {}
for key, value in coordination_numbers.items():
    for el in value:
        label_to_CN[el] = key

ligand_coords = {'lin': [[-1., 0, 0], [1, 0, 0]],
                 'trv': [[-0.5, -1*np.sqrt(3)/2, 0], [1, 0, 0]],
                 'tri': [[-0.5, -1*np.sqrt(3)/2, 0], [1, 0, 0], [-0.5, 1*np.sqrt(3)/2, 0]]}

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
trigonal_bipyramidal = [[0,0,1.], [0,0,-1]] + copy.copy(ligand_coords['tri'])

ligand_coords['tbp'] = copy.copy(trigonal_bipyramidal)
ligand_coords['bva'] = copy.copy(trigonal_bipyramidal[1:])
ligand_coords['bvp'] = copy.copy(trigonal_bipyramidal[:-1])

# square pyramidals
sq_py = [[1.,0.,0.], [-1.,0.,0.], [0.,1.,0.],[0.,-1.,0.],[0.,0.,1.]]
ligand_coords['spy'] = copy.copy(sq_py)
ligand_coords['pyv'] = copy.copy(sq_py[1:])

# tpv/tpr (trigonal prism)
def trigonal_prism():
    l = 1./np.sqrt(7)
    h = np.sqrt(12./7)

    v1 = [0, 2*l, 0.5*h]
    v2 = [-0.5*h, -1*l, 0.5*h]
    v3 = [0.5*h, -1*l, 0.5*h]
    v4 = [0, 2*l, -0.5*h]
    v5 = [-0.5*h, -1*l, -0.5*h]
    v6 = [0.5*h, -1*l, -0.5*h]

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

def check_coordination_geometry(ligands, center, geom):
    """
    Given coordinates, and geometry name, get rmsd.
    """
    coordination_number = len(ligands)
    coords = normalize_bond_lengths(ligands, center)
    #geom_coords = np.array(ligand_coords[geom] + [[0,0,0]], dtype=np.float_)
    #_, distances, _ = icp_no_translation(coords, geom_coords)
    #rmsd = np.sqrt(np.sum(distances**2)) / (coordination_number + 1)

    coords = coords[:-1,:]
    geom_coords = np.array(ligand_coords[geom], dtype=np.float_)
    rmsd, R, pts2 = best_match(coords, geom_coords)

    less_than_AT = rmsd < assignment_threshold[geom]
    less_than_DT = rmsd < distortion_threshold[geom]

    return rmsd, less_than_AT, less_than_DT

def find_coordination_geometry(ligands, center):
    """
    Given an array of coordinates for the ligands (N by 3), and given
    the coordinates of the center (1 by 3), determine the geometry that
    most closely matches it.
    """
    # get coordination number
    coordination_number = len(ligands)

    coords = normalize_bond_lengths(ligands, center)
    coords = coords[:-1,:]

    possibilities = []
    for geom in coordination_numbers[coordination_number]:
        geom_coords = np.array(ligand_coords[geom], dtype=np.float_)
        rmsd, _, _ = best_match(coords, geom_coords)
        less_than_AT = rmsd < assignment_threshold[geom]
        less_than_DT = rmsd < distortion_threshold[geom]
        if not less_than_AT:
            continue
        possibilities.append((geom, rmsd, less_than_AT, less_than_DT))

    if len(possibilities) == 0:
        print("WARNING: no possible assignments found.")
        return None, None, None

    min_rmsd = possibilities[0][1]
    min_geom = possibilities[0][0]
    distorted = not possibilities[0][-1]
    for geom, rmsd, _, l_DT in possibilities:
        if rmsd < min_rmsd:
            min_rmsd = rmsd
            min_geom = geom
            distorted = not l_DT

    if distorted is None:
        return None, None, None

    return min_geom, min_rmsd, distorted

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

def plot_geometry(coords, plot_zero=True):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.set_aspect('equal')
    for coord in coords:
        ax.scatter3D(coord[0], coord[1], coord[2], color='blue')

    if plot_zero:
        ax.scatter3D([0], [0], [0], color='red')
    plt.show()

def compare_two_geometries(geom_1, geom_2):
    """
    Given 3-letter names of two geometries, compare them and calculate RMSD.
    """
    geom_1_coords = np.array(ligand_coords[geom_1], dtype=np.float_)
    geom_2_coords = np.array(ligand_coords[geom_2], dtype=np.float_)

    rmsd, R, pts2 = best_match(geom_1_coords, geom_2_coords)
    return rmsd


def best_fit_rotation(pts1, pts2):
    """
    Given two sets of points, find the best fit rotation.

    Assumes translation is 0.
    """
    n_dim = pts1.shape[1]

    # rotation matrix
    H = np.dot(pts1.T, pts2)
    U, S, Vt = np.linalg.svd(H)
    R = np.dot(Vt.T, U.T)

    # special reflection case
    if np.linalg.det(R) < 0:
        Vt[n_dim-1,:] *= -1
        R = np.dot(Vt.T, U.T)

    return R

def distance(pts1, pts2):
    """
    Given two sets of points, return distance between each pair of points.
    """
    if pts1.shape != pts2.shape:
        print("WARNING: input point sets must have the same dimensions.")
        raise ValueError

    cn = len(pts1)
    d = []
    for i in range(cn):
        dist = np.linalg.norm(pts1[i,:] - pts2[i,:])
        d.append(dist)
    return np.array(d)

def best_match(pts1, pts2):
    """
    Given two sets of points, find the best matching between the two
    that gives the smallest rmsd.
    """
    if pts1.shape != pts2.shape:
        print("WARNING: input point sets must have the same dimensions.")
        raise ValueError

    cn = len(pts1)
    indices = list(range(cn))

    min_dist = np.inf
    min_matching = None
    best_rot = None
    for matching in permutations(indices, cn):
        idx = list(matching)
        pts2_p = pts2[idx,:]
        R = best_fit_rotation(pts1, pts2_p)
        pts1p = R.dot(pts1.T).T
        d = distance(pts1p, pts2_p)
        rmsd = np.sqrt(np.mean(d**2))
        if rmsd < min_dist:
            min_dist = rmsd
            min_matching = matching
            best_rot = R

    return min_dist, best_rot, pts2[list(min_matching),:]


# TODO: is there a better way to do this?
geom_thresholds = """lin 0.25881904510252074 0.25881904510252074
trv 0.25881904510252074 0.25881904510252074
tri 0.16910197872576282 0.2113248654051871
tev 0.16910197872576285 0.24732126143424196
spv 0.2113248654051871 0.24732126143424196
tet 0.14644660940672613 0.30290544652768614
spl 0.18301270189221933 0.3266756093699887
bva 0.14644660940672616 0.3266756093699887
bvp 0.09229595564125724 0.2166341981837778
pyv 0.09229595564125728 0.27059805007309845
spy 0.16369153687076268 0.20123643203933214
tbp 0.16369153687076268 0.18353003987052263
tpv 0.1835300398705227 0.20123643203933214
oct 0.1421883347398175 0.2966011848770199
tpr 0.16610333076733078 0.24934072472241658
pva 0.14306177797101235 0.2966011848770199
pvp 0.1421883347398175 0.24841972303098867
cof 0.14390499166611917 0.2023399366805026
con 0.17478572691393254 0.2787563034450018
ctf 0.1430617779710123 0.2152307878630878
ctn 0.14248638401127614 0.24934072472241653
pbp 0.12438818365078702 0.2751857287036503
coc 0.12719425753652283 0.2612199498875671
ctp 0.10847143014774503 0.256161519068391
hva 0.15655801081562157 0.2751857287036503
hvp 0.12438818365078703 0.24461665926091966
cuv 0.12719425753652286 0.188748023672228
sav 0.10847143014774507 0.2182230916888635
hbp 0.1464466094067262 0.25214083786810737
cub 0.11897933331668945 0.24842154942429004
sqa 0.10305711245039832 0.26560351903516544
boc 0.11897933331668949 0.21029799262993193
bts 0.10305711245039834 0.2567091783334269
btt 0.18115570978365955 0.26560351903516544
ttp 0.09716317741757498 0.09716317741757498
csa 0.09716317741757498 0.09716317741757498"""

distortion_threshold = {}
assignment_threshold = {}
for line in geom_thresholds.split('\n'):
    ls = line.split()
    g = ls[0]
    distortion_threshold[g] = float(ls[1])
    assignment_threshold[g] = float(ls[2])
