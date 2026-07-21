from gen.messages_pb2 import Molecule
from nodes.descriptors import descriptors
from gen.axiom_context import SecretStatus


class _TestContext:
    """Minimal AxiomContext implementation for unit tests."""

    class _Logger:
        def debug(self, msg: str, **attrs) -> None: pass
        def info(self, msg: str, **attrs) -> None: pass
        def warn(self, msg: str, **attrs) -> None: pass
        def error(self, msg: str, **attrs) -> None: pass

    class _Secrets:
        def __init__(self, m: dict, revoked: set) -> None:
            self._m = m or {}
            self._revoked = revoked or set()
        def get(self, name: str):
            v = self._m.get(name)
            return (v, True) if v is not None else ("", False)
        def status(self, name: str) -> SecretStatus:
            if name in self._m:
                return SecretStatus.AVAILABLE
            if name in self._revoked:
                return SecretStatus.REVOKED
            return SecretStatus.UNSET

    def __init__(self, secrets_map: dict | None = None, revoked_names: set | None = None) -> None:
        self.log = self._Logger()
        self.secrets = self._Secrets(secrets_map or {}, revoked_names)
        self.execution_id = "test-execution-id"
        self.flow_id = "test-flow-id"
        self.tenant_id = "test-tenant-id"


# Independent oracle: PubChem CID 2244 (aspirin) computed properties —
# MolecularFormula=C9H8O4, MolecularWeight=180.16, ExactMass=180.04225873,
# TPSA=63.6, HeavyAtomCount=13, HBondDonorCount=1 — fetched live from
# https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/property/...
# TPSA in particular is Ertl's fixed fragment-contribution algorithm, which
# both RDKit and PubChem's Cactvs toolkit implement identically, so an exact
# match is expected (not just "close").
def test_descriptors_aspirin_matches_pubchem():
    ax = _TestContext()
    result = descriptors(ax, Molecule(smiles="CC(=O)Oc1ccccc1C(=O)O"))
    assert result.valid is True
    assert result.canonical_smiles == "CC(=O)Oc1ccccc1C(=O)O"
    assert abs(result.molecular_weight - 180.16) < 0.01
    assert abs(result.exact_mass - 180.04225873) < 1e-4
    assert abs(result.tpsa - 63.6) < 0.01
    assert result.heavy_atom_count == 13
    assert result.hbd == 1
    assert result.ring_count == 1
    assert result.aromatic_ring_count == 1
    assert result.formal_charge == 0
    assert result.num_heteroatoms == 4  # 4 oxygens, no nitrogens


# hba/rotatable_bonds/logp/molar_refractivity are RDKit's own specific
# algorithm definitions (there is no single universal cross-tool standard
# for "H-bond acceptor" or "rotatable bond" the way there is for TPSA/mass;
# PubChem's Cactvs-based counts for aspirin are HBA=4, RotatableBondCount=3
# — genuinely different definitions, not a bug). These values are pinned
# against a direct, independent call to the underlying rdkit.Chem.Descriptors
# functions (not through this node's code path) for RDKit 2026.3.4.
def test_descriptors_aspirin_rdkit_specific_fields():
    from rdkit import Chem
    from rdkit.Chem import Descriptors as RDDescriptors

    mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
    expected_logp = RDDescriptors.MolLogP(mol)
    expected_hba = RDDescriptors.NumHAcceptors(mol)
    expected_rotb = RDDescriptors.NumRotatableBonds(mol)
    expected_mr = RDDescriptors.MolMR(mol)

    ax = _TestContext()
    result = descriptors(ax, Molecule(smiles="CC(=O)Oc1ccccc1C(=O)O"))
    assert result.logp == expected_logp
    assert result.hba == expected_hba
    assert result.rotatable_bonds == expected_rotb
    assert result.molar_refractivity == expected_mr


def test_descriptors_no_heteroatoms():
    # Ethane: pure C/H, no heteroatoms, no rings, no polar surface.
    ax = _TestContext()
    result = descriptors(ax, Molecule(smiles="CC"))
    assert result.valid is True
    assert result.num_heteroatoms == 0
    assert result.tpsa == 0.0
    assert result.ring_count == 0
    assert result.hbd == 0
    assert result.hba == 0


def test_descriptors_invalid_smiles_returns_structured_error():
    ax = _TestContext()
    result = descriptors(ax, Molecule(smiles="((("))
    assert result.valid is False
    assert result.error != ""
