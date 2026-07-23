# chem-tools

Composable Axiom nodes for **cheminformatics / molecular property computation**,
wrapping [RDKit](https://www.rdkit.org/) (BSD-3-Clause) — the de-facto standard
open-source cheminformatics toolkit that owns every algorithm here. Built for the
Axiom marketplace under the handle `christiangeorgelucas`.

Every node identifies a molecule by **SMILES** (or an MDL molblock) through the
shared `Molecule` envelope (or a lean node-specific message using the same field
vocabulary — `smiles` / `canonical_smiles` / `valid` / `error`). All nodes are
pure input→output: deterministic, no state, filesystem, network, or secrets.

## Use it from your agent or app

Every node in this package is a **live, auto-scaling API endpoint** on the
[Axiom](https://axiomide.com) marketplace — call it from an AI agent or your own
code, with nothing to self-host.

**📦 See it on the marketplace:**
https://dev.axiomide.com/marketplace/christiangeorgelucas/chem-tools@0.1.0

**Hook it up to an AI agent (MCP).** Add Axiom's hosted MCP server to any MCP
client and every node becomes a typed tool your agent can call — search the
catalog, inspect a schema, and invoke it directly.

```bash
# Claude Code
claude mcp add --transport http axiom https://api.axiomide.com/mcp \
  --header "Authorization: Bearer $AXIOM_API_KEY"
```

Claude Desktop, Cursor, or any config-based client:

```json
{
  "mcpServers": {
    "axiom": {
      "type": "http",
      "url": "https://api.axiomide.com/mcp",
      "headers": { "Authorization": "Bearer YOUR_AXIOM_API_KEY" }
    }
  }
}
```

**Call it from the CLI.**

```bash
axiom invoke christiangeorgelucas/chem-tools/Parse --input '{ ... }'
```

**Call it over HTTP.**

```bash
curl -X POST https://api.axiomide.com/invocations/v1/nodes/christiangeorgelucas/chem-tools/0.1.0/Parse \
  -H "Authorization: Bearer $AXIOM_API_KEY" \
  -H 'Content-Type: application/json' \
  -d '{ ... }'
```

> Input/output schema for each node is on the marketplace page above, or via
> `axiom inspect node christiangeorgelucas/chem-tools/Parse`.

### Get started free

Install the CLI:

```bash
# macOS / Linux — Homebrew
brew install axiomide/tap/axiom

# macOS / Linux — install script
curl -fsSL https://raw.githubusercontent.com/AxiomIDE/axiom-releases/main/install.sh | sh
```

**Windows:** download the `windows/amd64` `.zip` from the
[releases page](https://github.com/AxiomIDE/axiom-releases/releases), unzip it,
and put `axiom.exe` on your `PATH`.

Then `axiom version` to verify, `axiom login` (GitHub or Google) to authenticate,
and create an API key under **Console → API Keys**. Docs and sign-up at
**[axiomide.com](https://axiomide.com)**.

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
