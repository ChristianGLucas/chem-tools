from gen.messages_pb2 import FromInChIInput
from nodes.from_in_ch_i import from_in_ch_i
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


# Independent oracle: this is PubChem CID 2244's own published InChI string
# for aspirin (an independently-generated input, not produced by this
# package's ToInChI node) — decoding it must yield aspirin's known
# structure: formula C9H8O4, MW 180.16, 13 heavy atoms, 1 ring.
_PUBCHEM_ASPIRIN_INCHI = "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"


def test_from_inchi_pubchem_aspirin_string_decodes_correctly():
    ax = _TestContext()
    result = from_in_ch_i(ax, FromInChIInput(inchi=_PUBCHEM_ASPIRIN_INCHI))
    assert result.valid is True
    assert result.canonical_smiles == "CC(=O)Oc1ccccc1C(=O)O"
    assert result.formula == "C9H8O4"
    assert abs(result.molecular_weight - 180.16) < 0.01
    assert result.num_atoms == 13
    assert result.num_rings == 1


def test_from_inchi_round_trips_with_to_inchi():
    from gen.messages_pb2 import ToInChIInput
    from nodes.to_in_ch_i import to_in_ch_i

    ax = _TestContext()
    forward = to_in_ch_i(ax, ToInChIInput(smiles="c1ccc2ccccc2c1"))  # naphthalene
    backward = from_in_ch_i(ax, FromInChIInput(inchi=forward.inchi))
    assert backward.valid is True
    assert backward.canonical_smiles == "c1ccc2ccccc2c1"


def test_from_inchi_invalid_string_returns_structured_error():
    ax = _TestContext()
    result = from_in_ch_i(ax, FromInChIInput(inchi="not an inchi string"))
    assert result.valid is False
    assert result.error != ""


def test_from_inchi_empty_input_returns_structured_error():
    ax = _TestContext()
    result = from_in_ch_i(ax, FromInChIInput())
    assert result.valid is False
    assert result.error != ""
