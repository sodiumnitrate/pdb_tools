# pdb_tools
A Python package written to fetch protein structures from PDB and analyze the geometry of metal binding.
Inspired by [MetalPDB](https://metalpdb.cerm.unifi.it/).

Utilizes [GEMMI](https://github.com/project-gemmi/gemmi) to perform symmetry operations and other advanced geometry operations.
Uses [pypdb](https://github.com/williamgilpin/pypdb) to interface with PDB's GraphQL API.

## Installation
### From repository
```
git clone https://github.com/sodiumnitrate/pdb_tools.git
cd pdb_tools
pip install .
```

## Examples
Get a structure from PDB, select an atom, and find its neighboring atoms:
```
from pdb_tools.entry import Entry
from pdb_tools.utils import cra_from_atom

entry = Entry('6q6b')
atom = entry.select_atoms('(CU1)')[0]
neighboring_atoms = entry.find_neighboring_atoms(atom)
```

Find aromatic ring center and normal vector:
```
entry = Entry('1GST')
res = ('PHE', 'A', 56)
com, normal = entry.get_aromatic_ring_center_and_normal(res)
```

Find pairs of atoms within crystal contacts:
```
entry = Entry('6q6b')
results = entry.find_contact_pairs()
```

Find metal sites and assign coordination geometry based on neighbors:
```
entry = Entry('1a2v')
entry.get_metal_sites('Cu')

for site in entry.metal_sites:
    site.assign_geometry()
    print(f"geometry: {site.print_geometry()} (site.assigned_geometry), rmsd: {site.rmsd:.3f}")
    print("ligands:")
    for ligand in site.ligands:
        print(f"{ligand.residue.name}")
```

## Assigning coordination geometry
Loosely following [Andreini et al.](https://academic.oup.com/bioinformatics/article/28/12/1658/270149?login=false), the coordination geometry is assigned as follows:

1. Translate and scale coordinates such that metal center is at `(0,0,0)` and the distance between each ligand atom and metal center is `3A`.
2. Given the coordination number (number of ligands), compare with geometries with the same coordination number by finding the rotation of best fit.
3. Calculate the root mean squared distance (RMSD) from each candidate geometry and assign to the nearest geometry. Check geometry thresholds to make sure geometry can be assigned unambiguously.

## Limitations
Geometry finding is currently not optimized beyond coordination number `>=7`, meaning assignment takes too long or is intractable.
