"""Shared helpers for christiangeorgelucas/chem-tools nodes.

Centralizes SMILES/molblock parsing, the basic-facts extraction every node
that echoes a Molecule needs, and the structured-error contract every node
follows: a node never raises out to the caller — on malformed input it
returns its output message with computed fields empty, `valid=False`, and
`error` set to a short human-readable message.

"""
from __future__ import annotations

from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import Descriptors, rdMolDescriptors, rdFingerprintGenerator, MACCSkeys

# Silence RDKit's C++ stderr logging (parse warnings/errors) — we surface our
# own structured error messages instead of letting the native logger spam
# stderr on every malformed input.
rdBase.DisableLog("rdApp.*")

# Fingerprint bit-length bounds — 0 means "use default"; anything above the
# cap is rejected rather than silently clamped, so a caller who asked for
# something out of range gets a clear error instead of a quietly-different
# answer.
DEFAULT_MORGAN_RADIUS = 2
DEFAULT_N_BITS = 2048
MAX_RADIUS = 8
MIN_N_BITS = 16
MAX_N_BITS = 16_384

# GetSubstructMatches cap — bounds both the compute cost of enumerating
# matches on a highly symmetric molecule and the output size.
MAX_SUBSTRUCTURE_MATCHES = 100


class ChemToolsError(Exception):
    """Carries a short, human-readable error message for a node to surface
    verbatim in its output message's `error` field."""


def mol_from_input(smiles: str, molblock: str = "") -> "Chem.Mol":
    """Parse a molecule from a SMILES string (preferred) or an MDL molblock
    (used when smiles is empty). Raises ChemToolsError with a descriptive
    message on any failure — never returns None, never lets an RDKit
    exception escape uncaught.
    """
    if smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ChemToolsError(f"could not parse SMILES: {smiles!r}")
        return mol
    if molblock:
        mol = Chem.MolFromMolBlock(molblock)
        if mol is None:
            raise ChemToolsError("could not parse molblock")
        return mol
    raise ChemToolsError("either smiles or molblock is required")


def basic_facts(mol: "Chem.Mol") -> dict:
    """Return the core identity fields the Molecule envelope carries as
    output, as a plain kwargs dict — the caller constructs the message
    (`Molecule(**basic_facts(mol))`) rather than this helper doing it, so the
    node body is always the one building the typed output.
    """
    return {
        "canonical_smiles": Chem.MolToSmiles(mol),
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "molecular_weight": Descriptors.MolWt(mol),
        "num_atoms": mol.GetNumAtoms(),
        "num_bonds": mol.GetNumBonds(),
        "num_rings": rdMolDescriptors.CalcNumRings(mol),
        "valid": True,
        "error": "",
    }


def query_mol_from(query: str, query_is_smiles: bool) -> "Chem.Mol":
    if not query:
        raise ChemToolsError("query is required")
    if query_is_smiles:
        qmol = Chem.MolFromSmiles(query)
        if qmol is None:
            raise ChemToolsError(f"could not parse query SMILES: {query!r}")
        return qmol
    qmol = Chem.MolFromSmarts(query)
    if qmol is None:
        raise ChemToolsError(f"could not parse query SMARTS: {query!r}")
    return qmol


def resolve_fp_params(radius: int, n_bits: int) -> tuple[int, int]:
    """Apply defaults for radius/n_bits (0 => default) and enforce bounds.
    Raises ChemToolsError if an explicitly-requested value is out of range.
    """
    r = radius if radius else DEFAULT_MORGAN_RADIUS
    n = n_bits if n_bits else DEFAULT_N_BITS
    if r < 0 or r > MAX_RADIUS:
        raise ChemToolsError(f"radius must be between 0 and {MAX_RADIUS}, got {radius}")
    if n < MIN_N_BITS or n > MAX_N_BITS:
        raise ChemToolsError(
            f"n_bits must be between {MIN_N_BITS} and {MAX_N_BITS}, got {n_bits}"
        )
    return r, n


def compute_fingerprint(mol: "Chem.Mol", fp_type: int, radius: int, n_bits: int):
    """Return (ExplicitBitVect, resolved_n_bits) for the requested
    FingerprintType (proto enum int value: 0=MORGAN, 1=MACCS, 2=ATOM_PAIR,
    3=TOPOLOGICAL_TORSION).
    """
    if fp_type == 1:  # MACCS
        return MACCSkeys.GenMACCSKeys(mol), 167
    r, n = resolve_fp_params(radius, n_bits)
    if fp_type == 0:  # MORGAN
        gen = rdFingerprintGenerator.GetMorganGenerator(radius=r, fpSize=n)
    elif fp_type == 2:  # ATOM_PAIR
        gen = rdFingerprintGenerator.GetAtomPairGenerator(fpSize=n)
    elif fp_type == 3:  # TOPOLOGICAL_TORSION
        gen = rdFingerprintGenerator.GetTopologicalTorsionGenerator(fpSize=n)
    else:
        raise ChemToolsError(f"unknown fp_type: {fp_type}")
    return gen.GetFingerprint(mol), n


def fp_to_hex(fp) -> str:
    bitstring = fp.ToBitString()
    hexlen = (len(bitstring) + 3) // 4
    return format(int(bitstring, 2), f"0{hexlen}x")


def tanimoto(fp_a, fp_b) -> float:
    return DataStructs.TanimotoSimilarity(fp_a, fp_b)
