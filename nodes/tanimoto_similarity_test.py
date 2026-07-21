from gen.messages_pb2 import SimilarityInput
from nodes.tanimoto_similarity import tanimoto_similarity
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


# Independent oracle: Tanimoto(A, A) = |A n A| / |A u A| = 1.0 is a
# mathematical identity — true for ANY fingerprint algorithm, by definition
# of the Tanimoto/Jaccard coefficient, independent of RDKit's specific
# fingerprint implementation.
def test_tanimoto_self_similarity_is_exactly_one():
    ax = _TestContext()
    result = tanimoto_similarity(
        ax,
        SimilarityInput(
            smiles_a="CC(=O)Oc1ccccc1C(=O)O", smiles_b="CC(=O)Oc1ccccc1C(=O)O"
        ),
    )
    assert result.valid is True
    assert result.tanimoto == 1.0


# Independent oracle: different spellings of the identical molecule must
# also score exactly 1.0 (same structure -> same fingerprint).
def test_tanimoto_same_molecule_different_spelling_is_one():
    ax = _TestContext()
    result = tanimoto_similarity(
        ax, SimilarityInput(smiles_a="c1ccccc1", smiles_b="C1=CC=CC=C1")
    )
    assert result.valid is True
    assert result.tanimoto == 1.0


def test_tanimoto_dissimilar_molecules_below_one():
    ax = _TestContext()
    result = tanimoto_similarity(
        ax, SimilarityInput(smiles_a="C", smiles_b="CC(=O)Oc1ccccc1C(=O)O")
    )
    assert result.valid is True
    assert 0.0 <= result.tanimoto < 1.0


def test_tanimoto_is_symmetric():
    ax = _TestContext()
    r1 = tanimoto_similarity(
        ax, SimilarityInput(smiles_a="CCO", smiles_b="CC(=O)Oc1ccccc1C(=O)O")
    )
    r2 = tanimoto_similarity(
        ax, SimilarityInput(smiles_a="CC(=O)Oc1ccccc1C(=O)O", smiles_b="CCO")
    )
    assert r1.tanimoto == r2.tanimoto


def test_tanimoto_invalid_smiles_returns_structured_error():
    ax = _TestContext()
    result = tanimoto_similarity(
        ax, SimilarityInput(smiles_a="not_valid(((", smiles_b="CCO")
    )
    assert result.valid is False
    assert result.error != ""
