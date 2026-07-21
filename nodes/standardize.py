from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

from gen.messages_pb2 import StandardizeInput, StandardizeOutput
from gen.axiom_context import AxiomContext
from nodes._common import ChemToolsError, mol_from_input

_normalizer = rdMolStandardize.Normalizer()
_frag_chooser = rdMolStandardize.LargestFragmentChooser()
_uncharger = rdMolStandardize.Uncharger()


def standardize(ax: AxiomContext, input: StandardizeInput) -> StandardizeOutput:
    """Normalize a molecule for downstream comparison: keep only the largest
    organic fragment (strip counterions/solvates) and neutralize formal
    charges where a neutral protonation state exists. Reports whether
    standardization changed the structure.
    """
    try:
        mol = mol_from_input(input.smiles)
        original_canonical = Chem.MolToSmiles(mol)

        standardized = _normalizer.normalize(mol)
        if not input.keep_salts:
            standardized = _frag_chooser.choose(standardized)
        if not input.keep_charges:
            standardized = _uncharger.uncharge(standardized)

        result_smiles = Chem.MolToSmiles(standardized)
        return StandardizeOutput(
            smiles=result_smiles,
            changed=result_smiles != original_canonical,
            valid=True,
        )
    except ChemToolsError as e:
        return StandardizeOutput(valid=False, error=str(e))
    except Exception:  # noqa: BLE001 - malformed input must never crash the node
        return StandardizeOutput(valid=False, error="failed to standardize")
