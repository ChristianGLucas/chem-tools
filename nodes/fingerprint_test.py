from gen.messages_pb2 import FingerprintInput, FingerprintType
from nodes.fingerprint import fingerprint
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


# Independent oracle / structural invariant: a fingerprint is a function of
# the molecule's structure, not of how its SMILES was spelled — two
# different, hand-written SMILES for benzene must produce the IDENTICAL
# fingerprint. True by definition of a structural fingerprint, verifiable
# without trusting the implementation under test.
def test_fingerprint_structure_invariant_to_spelling():
    ax = _TestContext()
    r1 = fingerprint(ax, FingerprintInput(smiles="c1ccccc1"))
    r2 = fingerprint(ax, FingerprintInput(smiles="C1=CC=CC=C1"))
    assert r1.valid is True and r2.valid is True
    assert r1.fingerprint_hex == r2.fingerprint_hex


def test_fingerprint_morgan_default_params():
    ax = _TestContext()
    result = fingerprint(ax, FingerprintInput(smiles="CC(=O)Oc1ccccc1C(=O)O"))
    assert result.valid is True
    assert result.fp_type == FingerprintType.MORGAN
    assert result.n_bits == 2048
    assert result.bit_count > 0
    assert len(result.fingerprint_hex) == 2048 // 4


def test_fingerprint_maccs_is_fixed_167_bits():
    # MACCS is a fixed, published 167-bit key set (166 meaningful keys plus
    # an unused bit 0) — independent of any radius/n_bits request.
    ax = _TestContext()
    result = fingerprint(
        ax, FingerprintInput(smiles="CC(=O)Oc1ccccc1C(=O)O", fp_type=FingerprintType.MACCS, n_bits=999)
    )
    assert result.valid is True
    assert result.n_bits == 167


def test_fingerprint_different_molecules_differ():
    ax = _TestContext()
    a = fingerprint(ax, FingerprintInput(smiles="C"))  # methane
    b = fingerprint(ax, FingerprintInput(smiles="CC(=O)Oc1ccccc1C(=O)O"))  # aspirin
    assert a.fingerprint_hex != b.fingerprint_hex


def test_fingerprint_invalid_smiles_returns_structured_error():
    ax = _TestContext()
    result = fingerprint(ax, FingerprintInput(smiles="not_valid((("))
    assert result.valid is False
    assert result.error != ""


def test_fingerprint_out_of_range_n_bits_returns_structured_error():
    ax = _TestContext()
    result = fingerprint(ax, FingerprintInput(smiles="CCO", n_bits=999_999_999))
    assert result.valid is False
    assert result.error != ""
