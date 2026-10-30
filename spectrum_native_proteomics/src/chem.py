"""Monoisotopic chemistry constants and peptide mass utilities."""
import numpy as np

PROTON = 1.007276467
H2O = 18.010565
NH3 = 17.026549
# Residue monoisotopic masses; C is carbamidomethylated (fixed), 'm' = oxidised Met.
# I and L are isobaric and merged into L.
AA_MASS = {
    'G': 57.021464, 'A': 71.037114, 'S': 87.032028, 'P': 97.052764,
    'V': 99.068414, 'T': 101.047679, 'C': 160.030649, 'L': 113.084064,
    'N': 114.042927, 'D': 115.026943, 'Q': 128.058578, 'K': 128.094963,
    'E': 129.042593, 'M': 131.040485, 'H': 137.058912, 'F': 147.068414,
    'R': 156.101111, 'Y': 163.063329, 'W': 186.079313, 'm': 147.035400,
}
AAS = list(AA_MASS.keys())
# Approximate UniProt/Swiss-Prot amino-acid background frequencies (I folded into L)
AA_FREQ = {'A': .0825, 'R': .0553, 'N': .0406, 'D': .0545, 'C': .0137, 'Q': .0393,
           'E': .0675, 'G': .0707, 'H': .0227, 'L': .0965 + .0596, 'K': .0584,
           'M': .0242, 'F': .0386, 'P': .0470, 'S': .0656, 'T': .0534, 'W': .0108,
           'Y': .0292, 'V': .0687}

# token vocabulary: 0 PAD, 1 BOS, 2 EOS, 3.. amino acids
PAD, BOS, EOS = 0, 1, 2
TOK = {a: i + 3 for i, a in enumerate(AAS)}
ITOS = {v: k for k, v in TOK.items()}
VOCAB = len(AAS) + 3
TOK_MASS = np.zeros(VOCAB)
for a, i in TOK.items():
    TOK_MASS[i] = AA_MASS[a]


def peptide_mass(seq):
    """Neutral monoisotopic mass of a peptide."""
    return sum(AA_MASS[a] for a in seq) + H2O


def precursor_mz(seq, z):
    return (peptide_mass(seq) + z * PROTON) / z


def fragment_ladders(seq):
    """Singly protonated b and y ion m/z arrays (length n-1 each)."""
    r = np.array([AA_MASS[a] for a in seq])
    b = np.cumsum(r)[:-1] + PROTON
    y = np.cumsum(r[::-1])[:-1] + H2O + PROTON
    return b, y
