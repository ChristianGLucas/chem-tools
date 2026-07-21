from rdkit import Chem

from gen.messages_pb2 import CanonicalizeInput, CanonicalSmiles
from gen.axiom_context import AxiomContext
from nodes._common import ChemToolsError, mol_from_input


def canonicalize(ax: AxiomContext, input: CanonicalizeInput) -> CanonicalSmiles:
    """Reduce a SMILES string to RDKit's canonical form, with options to
    strip stereochemistry (non-isomeric canonical SMILES) or kekulize (write
    explicit single/double bonds instead of aromatic lowercase atoms). Two
    SMILES for the same molecule always canonicalize to the same string.
    """
    try:
        mol = mol_from_input(input.smiles)
        smiles = Chem.MolToSmiles(
            mol,
            isomericSmiles=not input.strip_stereo,
            kekuleSmiles=input.kekulize,
        )
        return CanonicalSmiles(canonical_smiles=smiles, valid=True)
    except ChemToolsError as e:
        return CanonicalSmiles(valid=False, error=str(e))
    except Exception:  # noqa: BLE001 - malformed input must never crash the node
        return CanonicalSmiles(valid=False, error="failed to canonicalize")
