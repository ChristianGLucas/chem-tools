from gen.messages_pb2 import SubstructureMatchInput
from nodes.substructure_match import substructure_match
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


# Independent oracle: aspirin (CC(=O)Oc1ccccc1C(=O)O) has TWO C(=O)O motifs
# by hand inspection of its structure — the acetate ester and the
# carboxylic acid — so a generic "C(=O)O" SMARTS query matches exactly
# twice. This is derivable from the structure by hand, not by running the
# code under test.
def test_substructure_match_carboxyl_like_group_in_aspirin():
    ax = _TestContext()
    result = substructure_match(
        ax, SubstructureMatchInput(smiles="CC(=O)Oc1ccccc1C(=O)O", query="C(=O)O")
    )
    assert result.valid is True
    assert result.matches is True
    assert result.match_count == 2
    assert len(result.match_details) == 2
    for m in result.match_details:
        assert len(m.atom_indices) == 3


# Negative control: aspirin has no sulfur at all, so a sulfonic-acid SMARTS
# cannot match — a hand-obvious chemistry fact.
def test_substructure_match_no_match():
    ax = _TestContext()
    result = substructure_match(
        ax,
        SubstructureMatchInput(smiles="CC(=O)Oc1ccccc1C(=O)O", query="S(=O)(=O)O"),
    )
    assert result.valid is True
    assert result.matches is False
    assert result.match_count == 0
    assert len(result.match_details) == 0


def test_substructure_match_query_is_smiles():
    ax = _TestContext()
    result = substructure_match(
        ax,
        SubstructureMatchInput(smiles="c1ccccc1O", query="c1ccccc1", query_is_smiles=True),
    )
    assert result.valid is True
    assert result.matches is True


def test_substructure_match_invalid_query_returns_structured_error():
    ax = _TestContext()
    result = substructure_match(
        ax, SubstructureMatchInput(smiles="CCO", query="not a valid smarts (((")
    )
    assert result.valid is False
    assert result.error != ""


def test_substructure_match_invalid_target_returns_structured_error():
    ax = _TestContext()
    result = substructure_match(ax, SubstructureMatchInput(smiles="not valid(((", query="C"))
    assert result.valid is False
    assert result.error != ""
