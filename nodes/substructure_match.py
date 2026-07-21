from gen.messages_pb2 import AtomMatch, SubstructureMatchInput, SubstructureMatchOutput
from gen.axiom_context import AxiomContext
from nodes._common import (
    MAX_SUBSTRUCTURE_MATCHES,
    ChemToolsError,
    mol_from_input,
    query_mol_from,
)


def substructure_match(ax: AxiomContext, input: SubstructureMatchInput) -> SubstructureMatchOutput:
    """Test whether a SMARTS (or SMILES) query pattern occurs within a
    target molecule. Returns whether it matches, the total match count, and
    up to 100 individual matches as target-atom-index sets in query-atom
    order.
    """
    try:
        mol = mol_from_input(input.smiles)
        query = query_mol_from(input.query, input.query_is_smiles)
        matches = mol.GetSubstructMatches(
            query,
            useChirality=input.use_chirality,
            maxMatches=MAX_SUBSTRUCTURE_MATCHES,
        )
        return SubstructureMatchOutput(
            matches=len(matches) > 0,
            match_count=len(matches),
            match_details=[AtomMatch(atom_indices=list(m)) for m in matches],
            valid=True,
        )
    except ChemToolsError as e:
        return SubstructureMatchOutput(valid=False, error=str(e))
    except Exception:  # noqa: BLE001 - malformed input must never crash the node
        return SubstructureMatchOutput(valid=False, error="failed to match substructure")
