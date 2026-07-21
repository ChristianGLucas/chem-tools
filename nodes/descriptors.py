from rdkit import Chem
from rdkit.Chem import Descriptors as RDDescriptors
from rdkit.Chem import rdMolDescriptors

from gen.messages_pb2 import Molecule, Descriptors
from gen.axiom_context import AxiomContext
from nodes._common import ChemToolsError, mol_from_input


def descriptors(ax: AxiomContext, input: Molecule) -> Descriptors:
    """Compute a standard batch of RDKit molecular descriptors for
    drug-likeness / QSAR use: molecular weight, exact mass, Crippen LogP,
    TPSA, H-bond donor/acceptor counts, rotatable bonds, ring and
    aromatic-ring counts, heavy-atom count, fraction Csp3, molar
    refractivity, formal charge, and heteroatom count.
    """
    try:
        mol = mol_from_input(input.smiles, input.molblock)
        num_heteroatoms = sum(
            1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in (1, 6)
        )
        return Descriptors(
            smiles=input.smiles,
            canonical_smiles=Chem.MolToSmiles(mol),
            molecular_weight=RDDescriptors.MolWt(mol),
            exact_mass=RDDescriptors.ExactMolWt(mol),
            logp=RDDescriptors.MolLogP(mol),
            tpsa=RDDescriptors.TPSA(mol),
            hbd=RDDescriptors.NumHDonors(mol),
            hba=RDDescriptors.NumHAcceptors(mol),
            rotatable_bonds=RDDescriptors.NumRotatableBonds(mol),
            ring_count=rdMolDescriptors.CalcNumRings(mol),
            aromatic_ring_count=rdMolDescriptors.CalcNumAromaticRings(mol),
            heavy_atom_count=mol.GetNumHeavyAtoms(),
            fraction_csp3=rdMolDescriptors.CalcFractionCSP3(mol),
            molar_refractivity=RDDescriptors.MolMR(mol),
            formal_charge=Chem.GetFormalCharge(mol),
            num_heteroatoms=num_heteroatoms,
            valid=True,
        )
    except ChemToolsError as e:
        return Descriptors(valid=False, error=str(e))
    except Exception:  # noqa: BLE001 - malformed input must never crash the node
        return Descriptors(valid=False, error="failed to compute descriptors")
