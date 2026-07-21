from gen.messages_pb2 import Molecule
from gen.axiom_context import AxiomContext
from nodes._common import ChemToolsError, basic_facts, mol_from_input


def parse(ax: AxiomContext, input: Molecule) -> Molecule:
    """Parse a SMILES string (or MDL molblock, when smiles is empty) into a
    validated molecule: canonical isomeric SMILES, Hill-notation molecular
    formula, average molecular weight, and heavy-atom/bond/ring counts.
    Malformed input returns valid=false with a structured error instead of
    raising.
    """
    try:
        mol = mol_from_input(input.smiles, input.molblock)
        return Molecule(**basic_facts(mol))
    except ChemToolsError as e:
        return Molecule(valid=False, error=str(e))
    except Exception:  # noqa: BLE001 - malformed input must never crash the node
        return Molecule(valid=False, error="failed to parse molecule")
