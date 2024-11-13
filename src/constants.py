from collections import defaultdict

APP_VERSION = 1.0

ERROR_SIM_NOT_LOADED_MSG = "Simulation hasn't been loaded (or load failed)"
ERR_SEGMENT_NOT_FOUND = "Segment '{}' is not present in the simulation/results"

# In confidently determining atomic element from its mass:
# abs(elem_to_determine.mass - elem_possible.mass) must be <= ATOMTABLE_ELEMENT_MASS_CONFIDENCE_THRESHOLD
ATOMTABLE_ELEMENT_MASS_CONFIDENCE_THRESHOLD = 0.01

# Masses of atoms given their element
# from: https://gist.github.com/lukasrichters14/c862644d4cbcf2d67252a484b7c6049c
ATOMTABLE_ELEMENTS_MASSES = {'H' : 1.008, 'He' : 4.003, 'Li' : 6.941, 'Be' : 9.012,
                             'B' : 10.811, 'C' : 12.011, 'N' : 14.007, 'O' : 15.999,
                             'F' : 18.998, 'Ne' : 20.180, 'Na' : 22.990, 'Mg' : 24.305,
                             'Al' : 26.982, 'Si' : 28.086, 'P' : 30.974, 'S' : 32.066,
                             'Cl' : 35.453, 'Ar' : 39.948, 'K' : 39.098, 'Ca' : 40.078,
                             'Sc' : 44.956, 'Ti' : 47.867, 'V' : 50.942, 'Cr' : 51.996,
                             'Mn' : 54.938, 'Fe' : 55.845, 'Co' : 58.933, 'Ni' : 58.693,
                             'Cu' : 63.546, 'Zn' : 65.38, 'Ga' : 69.723, 'Ge' : 72.631,
                             'As' : 74.922, 'Se' : 78.971, 'Br' : 79.904, 'Kr' : 84.798,
                             'Rb' : 84.468, 'Sr' : 87.62, 'Y' : 88.906, 'Zr' : 91.224,
                             'Nb' : 92.906, 'Mo' : 95.95, 'Tc' : 98.907, 'Ru' : 101.07,
                             'Rh' : 102.906, 'Pd' : 106.42, 'Ag' : 107.868, 'Cd' : 112.414,
                             'In' : 114.818, 'Sn' : 118.711, 'Sb' : 121.760, 'Te' : 126.7,
                             'I' : 126.904, 'Xe' : 131.294, 'Cs' : 132.905, 'Ba' : 137.328,
                             'La' : 138.905, 'Ce' : 140.116, 'Pr' : 140.908, 'Nd' : 144.243,
                             'Pm' : 144.913, 'Sm' : 150.36, 'Eu' : 151.964, 'Gd' : 157.25,
                             'Tb' : 158.925, 'Dy': 162.500, 'Ho' : 164.930, 'Er' : 167.259,
                             'Tm' : 168.934, 'Yb' : 173.055, 'Lu' : 174.967, 'Hf' : 178.49,
                             'Ta' : 180.948, 'W' : 183.84, 'Re' : 186.207, 'Os' : 190.23,
                             'Ir' : 192.217, 'Pt' : 195.085, 'Au' : 196.967, 'Hg' : 200.592,
                             'Tl' : 204.383, 'Pb' : 207.2, 'Bi' : 208.980, 'Po' : 208.982,
                             'At' : 209.987, 'Rn' : 222.081, 'Fr' : 223.020, 'Ra' : 226.025,
                             'Ac' : 227.028, 'Th' : 232.038, 'Pa' : 231.036, 'U' : 238.029,
                             'Np' : 237, 'Pu' : 244, 'Am' : 243, 'Cm' : 247, 'Bk' : 247,
                             'Ct' : 251, 'Es' : 252, 'Fm' : 257, 'Md' : 258, 'No' : 259,
                             'Lr' : 262, 'Rf' : 261, 'Db' : 262, 'Sg' : 266, 'Bh' : 264,
                             'Hs' : 269, 'Mt' : 268, 'Ds' : 271, 'Rg' : 272, 'Cn' : 285,
                             'Nh' : 284, 'Fl' : 289, 'Mc' : 288, 'Lv' : 292, 'Ts' : 294,
                             'Og' : 294}

# number of expected bonds for chemical element
# tuple of (neutral, [possible configurations])
_EXPECTED_BONDS_SRC: dict[str, tuple[int, list[int]]] = {
    'C': (4, [4]),
    'H': (1, [1]),
    'O': (2, [1, 2]),
    'N': (3, [2, 3, 4]),
    'P': (5, [5]),
    'S': (2, [1, 2, 3]),
    'F': (1, [1]),
    'Cl': (1, [1]),
    'Br': (1, [1]),
    'I': (1, [1]),
    'Mg': (1, [1, 2]),
    'Au': (0, [1, 3])
}
ATOMTABLE_ELEMENT_EXPECTED_BONDS_NEUTRAL = 0
ATOMTABLE_ELEMENT_EXPECTED_BONDS_POSSIBLE = 1

ATOMTABLE_ELEMENT_EXPECTED_BONDS = defaultdict(lambda: (0, [0]))
for elem, bond_info in _EXPECTED_BONDS_SRC.items():
    ATOMTABLE_ELEMENT_EXPECTED_BONDS[elem] = bond_info

# if approximating atomic element from its mass, number of closest results
RESOLVER_ATOM_APPROX_FIND_SIMILAR_N_MASSES = 5

# all recognized amino acids' 3-letter codes and their corresponding 1-letter code
KNOWN_AMINO_ACIDS: dict[str, str] = {
    'ALA': 'A',
    'ARG': 'R',
    'ASN': 'N',
    'ASP': 'D',
    'CYS': 'C',
    'GLU': 'E',
    'GLN': 'Q',
    'GLY': 'G',
    'HIS': 'H',
    'HYP': 'O',
    'ILE': 'I',
    'LEU': 'L',
    'LYS': 'K',
    'MET': 'M',
    'PHE': 'F',
    'PRO': 'P',
    'SER': 'S',
    'THR': 'T',
    'TRP': 'W',
    'TYR': 'Y',
    'VAL': 'V', }


NETWORKX_COLORMAP_ATOMS = defaultdict(lambda: "gray")
for element, color in {
    'H': 'white',
    'C': 'black',
    'N': 'green',
    'O': 'blue',
    'P': 'red',
    'S': 'yellow',
}.items():
    NETWORKX_COLORMAP_ATOMS[element] = color