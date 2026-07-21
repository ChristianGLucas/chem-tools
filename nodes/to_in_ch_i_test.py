from gen.messages_pb2 import ToInChIInput
from nodes.to_in_ch_i import to_in_ch_i
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


# Independent oracle: PubChem CID 2244 (aspirin) publishes both InChI and
# InChIKey — fetched live from
# https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/property/InChI,InChIKey/JSON
#   InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)
#   InChIKey: BSYNRYMUTXBXSQ-UHFFFAOYSA-N
# InChI is a fixed IUPAC-standard algorithm, so an exact string match is
# the correct bar (not "close enough").
def test_to_inchi_aspirin_matches_pubchem():
    ax = _TestContext()
    result = to_in_ch_i(ax, ToInChIInput(smiles="CC(=O)Oc1ccccc1C(=O)O"))
    assert result.valid is True
    assert result.inchi == "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"
    assert result.inchi_key == "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"


def test_to_inchi_invalid_smiles_returns_structured_error():
    ax = _TestContext()
    result = to_in_ch_i(ax, ToInChIInput(smiles="not_valid((("))
    assert result.valid is False
    assert result.error != ""


def test_to_inchi_is_deterministic():
    ax = _TestContext()
    r1 = to_in_ch_i(ax, ToInChIInput(smiles="CCO"))
    r2 = to_in_ch_i(ax, ToInChIInput(smiles="CCO"))
    assert r1.inchi == r2.inchi
    assert r1.inchi_key == r2.inchi_key
