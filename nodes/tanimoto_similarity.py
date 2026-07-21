from gen.messages_pb2 import SimilarityInput, SimilarityOutput
from gen.axiom_context import AxiomContext
from nodes._common import ChemToolsError, compute_fingerprint, mol_from_input, tanimoto


def tanimoto_similarity(ax: AxiomContext, input: SimilarityInput) -> SimilarityOutput:
    """Compute the Tanimoto (Jaccard) similarity coefficient between two
    molecules' fingerprints, in [0, 1]. The same fingerprint algorithm and
    parameters are applied to both molecules before comparing.
    """
    try:
        mol_a = mol_from_input(input.smiles_a)
        mol_b = mol_from_input(input.smiles_b)
        fp_a, _ = compute_fingerprint(mol_a, input.fp_type, input.radius, input.n_bits)
        fp_b, _ = compute_fingerprint(mol_b, input.fp_type, input.radius, input.n_bits)
        return SimilarityOutput(tanimoto=tanimoto(fp_a, fp_b), valid=True)
    except ChemToolsError as e:
        return SimilarityOutput(valid=False, error=str(e))
    except Exception:  # noqa: BLE001 - malformed input must never crash the node
        return SimilarityOutput(valid=False, error="failed to compute similarity")
