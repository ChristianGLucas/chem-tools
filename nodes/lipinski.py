from rdkit import Chem
from rdkit.Chem import Descriptors as RDDescriptors

from gen.messages_pb2 import Molecule, LipinskiOutput
from gen.axiom_context import AxiomContext
from nodes._common import ChemToolsError, mol_from_input


def lipinski(ax: AxiomContext, input: Molecule) -> LipinskiOutput:
    """Screen a molecule against Lipinski's Rule of Five (drug-likeness):
    molecular weight <= 500, LogP <= 5, H-bond donors <= 5, H-bond
    acceptors <= 10. Reports the per-rule pass/fail, violation count, and
    overall pass (<= 1 violation, the standard allowance).
    """
    try:
        mol = mol_from_input(input.smiles, input.molblock)
        mw = RDDescriptors.MolWt(mol)
        logp = RDDescriptors.MolLogP(mol)
        hbd = RDDescriptors.NumHDonors(mol)
        hba = RDDescriptors.NumHAcceptors(mol)

        mw_pass = mw <= 500
        logp_pass = logp <= 5
        hbd_pass = hbd <= 5
        hba_pass = hba <= 10
        violations = sum(1 for p in (mw_pass, logp_pass, hbd_pass, hba_pass) if not p)

        return LipinskiOutput(
            canonical_smiles=Chem.MolToSmiles(mol),
            molecular_weight=mw,
            logp=logp,
            hbd=hbd,
            hba=hba,
            mw_pass=mw_pass,
            logp_pass=logp_pass,
            hbd_pass=hbd_pass,
            hba_pass=hba_pass,
            violations=violations,
            passes=violations <= 1,
            valid=True,
        )
    except ChemToolsError as e:
        return LipinskiOutput(valid=False, error=str(e))
    except Exception:  # noqa: BLE001 - malformed input must never crash the node
        return LipinskiOutput(valid=False, error="failed to screen molecule")
