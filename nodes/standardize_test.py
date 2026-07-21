from gen.messages_pb2 import StandardizeInput
from nodes.standardize import standardize
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


# Independent oracle: sodium acetate (CC(=O)O.[Na+]) is a basic-chemistry
# textbook fact — stripping the sodium counterion and neutralizing leaves
# plain acetic acid, CC(=O)O. Hand-verifiable without RDKit.
def test_standardize_strips_salt_and_neutralizes():
    ax = _TestContext()
    result = standardize(ax, StandardizeInput(smiles="CC(=O)O.[Na+]"))
    assert result.valid is True
    assert result.smiles == "CC(=O)O"
    assert result.changed is True


def test_standardize_keep_salts_preserves_fragment():
    ax = _TestContext()
    result = standardize(ax, StandardizeInput(smiles="CC(=O)O.[Na+]", keep_salts=True))
    assert result.valid is True
    assert "." in result.smiles  # multi-component: fragment retained


def test_standardize_already_clean_molecule_is_unchanged():
    # Ethanol: already neutral, already a single fragment — nothing to do.
    ax = _TestContext()
    result = standardize(ax, StandardizeInput(smiles="CCO"))
    assert result.valid is True
    assert result.smiles == "CCO"
    assert result.changed is False


def test_standardize_invalid_smiles_returns_structured_error():
    ax = _TestContext()
    result = standardize(ax, StandardizeInput(smiles="not_valid((("))
    assert result.valid is False
    assert result.error != ""
