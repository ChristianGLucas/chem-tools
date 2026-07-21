# chem-tools

Composable Axiom nodes for **cheminformatics / molecular property computation**,
wrapping [RDKit](https://www.rdkit.org/) (BSD-3-Clause) — the de-facto standard
open-source cheminformatics toolkit that owns every algorithm here. Built for the
Axiom marketplace under the handle `christiangeorgelucas`.

Every node identifies a molecule by **SMILES** (or an MDL molblock) through the
shared `Molecule` envelope (or a lean node-specific message using the same field
vocabulary — `smiles` / `canonical_smiles` / `valid` / `error`). All nodes are
pure input→output: deterministic, no state, filesystem, network, or secrets.

## Nodes

| Node | Does |
|---|---|
| **Parse** | Validate a SMILES/molblock: canonical SMILES, formula, MW, atom/bond/ring counts. |
| **Descriptors** | Standard RDKit descriptor batch: MW, exact mass, LogP, TPSA, HBD/HBA, rotatable bonds, ring counts, fraction Csp3, molar refractivity, formal charge, heteroatom count. |
| **Canonicalize** | Reduce a SMILES to RDKit's canonical form (optionally strip stereo / kekulize). |
| **Standardize** | Strip salts/solvates to the largest fragment and neutralize charges. |
| **SubstructureMatch** | Test a SMARTS/SMILES query against a target molecule; returns match count + atom indices. |
| **Fingerprint** | Morgan/ECFP, MACCS, Atom-Pair, or Topological-Torsion fingerprint as hex. |
| **TanimotoSimilarity** | Tanimoto coefficient between two molecules' fingerprints. |
| **Lipinski** | Rule-of-Five drug-likeness screen with per-rule pass/fail. |
| **ToInChI** | Molecule → standard InChI string + InChIKey. |
| **FromInChI** | InChI string → molecule (canonical SMILES, formula, MW, counts). |
| **MurckoScaffold** | Extract the Murcko ring scaffold and its generic (all-carbon) framework. |

## Safety & determinism

- **Bounded input.** SMILES/molblock/SMARTS/InChI strings are length-capped
  before any RDKit call; fingerprint bit-length and Morgan radius are bounded;
  substructure matches are capped at 100 results — all before allocation, not
  after.
- **Structured errors.** Malformed input never raises — every node returns its
  output message with `valid=false` and a human-readable `error`.
- **No drawing/rendering nodes.** RDKit's `Draw` module bundles fonts under an
  ambiguous "found on the internet, no explicit license" notice; this package
  sticks to RDKit's core cheminformatics algorithms (parsing, descriptors,
  fingerprints, InChI) to keep the license tree unambiguous.

## License

MIT (this package). Wraps RDKit, BSD-3-Clause.
