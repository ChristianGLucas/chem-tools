from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold as RDMurckoScaffold

from gen.messages_pb2 import Molecule, ScaffoldOutput
from gen.axiom_context import AxiomContext
from nodes._common import ChemToolsError, mol_from_input


def murcko_scaffold(ax: AxiomContext, input: Molecule) -> ScaffoldOutput:
    """Extract the Murcko scaffold of a molecule (rings + connecting
    linkers, side chains removed) and its generic framework (scaffold with
    every atom set to carbon and every bond set to single) — used for
    scaffold-based clustering and analysis in compound libraries.
    """
    try:
        mol = mol_from_input(input.smiles, input.molblock)
        scaffold = RDMurckoScaffold.GetScaffoldForMol(mol)
        generic = RDMurckoScaffold.MakeScaffoldGeneric(scaffold)
        return ScaffoldOutput(
            canonical_smiles=Chem.MolToSmiles(mol),
            scaffold_smiles=Chem.MolToSmiles(scaffold),
            generic_scaffold_smiles=Chem.MolToSmiles(generic),
            valid=True,
        )
    except ChemToolsError as e:
        return ScaffoldOutput(valid=False, error=str(e))
    except Exception:  # noqa: BLE001 - malformed input must never crash the node
        return ScaffoldOutput(valid=False, error="failed to extract scaffold")
