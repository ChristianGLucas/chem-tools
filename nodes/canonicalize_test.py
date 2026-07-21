from gen.messages_pb2 import CanonicalizeInput
from nodes.canonicalize import canonicalize
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


# Independent oracle: canonicalization must be a function of the MOLECULE,
# not of the input string's spelling — two different, hand-written SMILES
# for benzene (aromatic lowercase form vs explicit Kekule form) must
# canonicalize to the identical string. This is a structural invariant true
# by the definition of "canonical form", not something derived from this
# package's own code.
def test_canonicalize_same_molecule_different_spellings_converge():
    ax = _TestContext()
    r1 = canonicalize(ax, CanonicalizeInput(smiles="c1ccccc1"))
    r2 = canonicalize(ax, CanonicalizeInput(smiles="C1=CC=CC=C1"))
    assert r1.valid is True and r2.valid is True
    assert r1.canonical_smiles == r2.canonical_smiles == "c1ccccc1"


def test_canonicalize_strip_stereo():
    # L-alanine has one defined stereocenter (@@). Default preserves it;
    # strip_stereo removes the marker but keeps the same connectivity.
    ax = _TestContext()
    with_stereo = canonicalize(ax, CanonicalizeInput(smiles="C[C@@H](N)C(=O)O"))
    without_stereo = canonicalize(
        ax, CanonicalizeInput(smiles="C[C@@H](N)C(=O)O", strip_stereo=True)
    )
    assert "@" in with_stereo.canonical_smiles
    assert "@" not in without_stereo.canonical_smiles
    assert without_stereo.canonical_smiles == "CC(N)C(=O)O"


def test_canonicalize_kekulize():
    ax = _TestContext()
    aromatic = canonicalize(ax, CanonicalizeInput(smiles="c1ccccc1"))
    kekulized = canonicalize(ax, CanonicalizeInput(smiles="c1ccccc1", kekulize=True))
    assert aromatic.canonical_smiles == "c1ccccc1"
    assert kekulized.canonical_smiles == "C1=CC=CC=C1"


def test_canonicalize_invalid_smiles_returns_structured_error():
    ax = _TestContext()
    result = canonicalize(ax, CanonicalizeInput(smiles="not_valid((("))
    assert result.valid is False
    assert result.error != ""
