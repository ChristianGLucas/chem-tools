from gen.messages_pb2 import FingerprintInput, FingerprintOutput
from gen.axiom_context import AxiomContext
from nodes._common import ChemToolsError, compute_fingerprint, fp_to_hex, mol_from_input
from rdkit import Chem


def fingerprint(ax: AxiomContext, input: FingerprintInput) -> FingerprintOutput:
    """Compute a molecular fingerprint — Morgan/ECFP (radius + folded bit
    length configurable), MACCS keys (fixed 167 bits), Atom-Pair, or
    Topological-Torsion — as a hex-encoded bit vector plus its popcount.
    Used for similarity search and clustering.
    """
    try:
        mol = mol_from_input(input.smiles)
        fp, n_bits = compute_fingerprint(mol, input.fp_type, input.radius, input.n_bits)
        return FingerprintOutput(
            canonical_smiles=Chem.MolToSmiles(mol),
            fp_type=input.fp_type,
            fingerprint_hex=fp_to_hex(fp),
            n_bits=n_bits,
            bit_count=fp.GetNumOnBits(),
            valid=True,
        )
    except ChemToolsError as e:
        return FingerprintOutput(valid=False, error=str(e))
    except Exception:  # noqa: BLE001 - malformed input must never crash the node
        return FingerprintOutput(valid=False, error="failed to compute fingerprint")
