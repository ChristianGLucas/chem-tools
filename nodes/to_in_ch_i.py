from rdkit.Chem import inchi

from gen.messages_pb2 import ToInChIInput, InChIOutput
from gen.axiom_context import AxiomContext
from nodes._common import ChemToolsError, mol_from_input


def to_in_ch_i(ax: AxiomContext, input: ToInChIInput) -> InChIOutput:
    """Convert a molecule to its standard InChI string and InChIKey (the
    fixed-length hashed identifier), via IUPAC's InChI algorithm as bundled
    in RDKit.
    """
    try:
        mol = mol_from_input(input.smiles)
        inchi_str = inchi.MolToInchi(mol, options=input.options)
        if not inchi_str:
            raise ChemToolsError("RDKit could not generate an InChI for this molecule")
        inchi_key = inchi.InchiToInchiKey(inchi_str)
        return InChIOutput(inchi=inchi_str, inchi_key=inchi_key, valid=True)
    except ChemToolsError as e:
        return InChIOutput(valid=False, error=str(e))
    except Exception:  # noqa: BLE001 - malformed input must never crash the node
        return InChIOutput(valid=False, error="failed to generate InChI")
