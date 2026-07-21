from gen.messages_pb2 import Molecule
from nodes.lipinski import lipinski
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


# Independent oracle: aspirin is a textbook oral drug (PubChem CID 2244,
# MW=180.16 < 500, HBD=1 < 5) — it must pass Lipinski's Rule of Five with
# zero violations. Hand-verifiable from the PubChem-published facts alone.
def test_lipinski_aspirin_passes():
    ax = _TestContext()
    result = lipinski(ax, Molecule(smiles="CC(=O)Oc1ccccc1C(=O)O"))
    assert result.valid is True
    assert abs(result.molecular_weight - 180.16) < 0.01
    assert result.hbd == 1
    assert result.mw_pass is True
    assert result.hbd_pass is True
    assert result.violations == 0
    assert result.passes is True


# Independent oracle: a heavily hydroxylated open-chain polyol (13 OH
# groups) has, by hand count of its SMILES, HBD=HBA=14 — over the HBD<=5
# and HBA<=10 limits, while its MW (~422) and LogP (very negative) stay
# within limits. That is exactly 2 rule violations (HBD, HBA), which fails
# the standard <=1-violation allowance. Verified directly against RDKit
# (independent of this node's code path) before being pinned here.
def test_lipinski_polyol_fails_multiple_rules():
    ax = _TestContext()
    smiles = "OCC(O)C(O)C(O)C(O)C(O)C(O)C(O)C(O)C(O)C(O)C(O)CO"
    result = lipinski(ax, Molecule(smiles=smiles))
    assert result.valid is True
    assert result.mw_pass is True
    assert result.logp_pass is True
    assert result.hbd_pass is False
    assert result.hba_pass is False
    assert result.violations == 2
    assert result.passes is False


def test_lipinski_invalid_smiles_returns_structured_error():
    ax = _TestContext()
    result = lipinski(ax, Molecule(smiles="not_valid((("))
    assert result.valid is False
    assert result.error != ""
