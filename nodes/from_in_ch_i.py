from rdkit.Chem import inchi

from gen.messages_pb2 import FromInChIInput, Molecule
from gen.axiom_context import AxiomContext
from nodes._common import ChemToolsError, basic_facts


def from_in_ch_i(ax: AxiomContext, input: FromInChIInput) -> Molecule:
    """Parse an InChI string back into a molecule: canonical SMILES,
    formula, molecular weight, and atom/bond/ring counts — the inverse of
    ToInChI, useful for round-tripping and InChI-keyed lookups.
    """
    try:
        if not input.inchi:
            raise ChemToolsError("inchi is required")
        mol = inchi.MolFromInchi(input.inchi)
        if mol is None:
            raise ChemToolsError(f"could not parse InChI: {input.inchi!r}")
        return Molecule(**basic_facts(mol))
    except ChemToolsError as e:
        return Molecule(valid=False, error=str(e))
    except Exception:  # noqa: BLE001 - malformed input must never crash the node
        return Molecule(valid=False, error="failed to parse InChI")
