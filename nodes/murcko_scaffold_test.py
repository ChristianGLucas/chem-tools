from gen.messages_pb2 import Molecule
from nodes.murcko_scaffold import murcko_scaffold
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


# Independent oracle: by the textbook definition of a Murcko scaffold
# (rings + linkers only, side chains/substituents removed), aspirin's
# scaffold is nothing but its single benzene ring — the acetoxy and
# carboxyl groups are terminal substituents, not linkers between rings.
# Hand-derivable from the structure, not from running RDKit.
def test_murcko_scaffold_aspirin_is_benzene():
    ax = _TestContext()
    result = murcko_scaffold(ax, Molecule(smiles="CC(=O)Oc1ccccc1C(=O)O"))
    assert result.valid is True
    assert result.scaffold_smiles == "c1ccccc1"
    assert result.generic_scaffold_smiles == "C1CCCCC1"


def test_murcko_scaffold_no_rings_is_empty():
    # A molecule with no ring system has an empty Murcko scaffold.
    ax = _TestContext()
    result = murcko_scaffold(ax, Molecule(smiles="CCCCCC"))
    assert result.valid is True
    assert result.scaffold_smiles == ""


def test_murcko_scaffold_invalid_smiles_returns_structured_error():
    ax = _TestContext()
    result = murcko_scaffold(ax, Molecule(smiles="not_valid((("))
    assert result.valid is False
    assert result.error != ""
