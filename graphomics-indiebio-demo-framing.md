# Graphomics — IndieBio Demo Framing & Feasibility Hunt

**Audience:** Vinith, Varun, Demetrius
**Occasion:** IndieBio demo, onsite SF, June 9 (held earlier per scheduling)
**Status:** Framing locked. Feasibility unverified — open questions at the end need answers before we commit to a *live* vs *recorded* demo.

-----

## 1. The thesis (the “why”)

Every pharma has a lossy, expensive translation gap between the people who **allocate capital** (CEO / portfolio) and the people who **defend the science** (CSO / bench). These two sit on opposite sides of the drug-development pipeline: the CEO lives on the right (portfolio, Go/No-Go gates, filings, ROI), the CSO’s truth lives on the left (target validity, mechanism, assay integrity). Every Go decision is a moment where capital meets science — and right now the translation is lossy in both directions. The CEO can’t see why a target is weak; the CSO can’t always express scientific risk in the language of portfolio consequence.

**Graphomics is the mechanistic-science substrate that closes that gap.** A federated, tiered, provenance-carrying knowledge layer that lets *any* agentic harness reason from a CEO-altitude question *down* to root-level mechanism and back *up* — without flattening the science on the way. Someone has to pay for the science, and the science has to lead to drugs that make that money back. We are the layer that keeps those two facts in one conversation.

-----

## 2. Competitive frame (Perceptic)

Perceptic is **not** a head-to-head competitor; it operates in the opposite direction.

|                      |Perceptic                                   |Graphomics                                        |
|----------------------|--------------------------------------------|--------------------------------------------------|
|Direction of reasoning|Top-down (reasons *down* from the portfolio)|Bottom-up (reasons *up* from the mechanism)       |
|Altitude served       |CEO / program decisions                     |The connective tissue across CEO↔CSO              |
|Treatment of science  |Inherits it as ground truth                 |Grounds it, scores it, makes it traceable         |
|Anchored in           |Existing org data + decisions               |Root-level mechanistic biology (MapForge/MicroMap)|
|Atlas module scope    |Clinical trial data                         |Everything upstream of the clinic                 |

One-liner if it comes up: *“Perceptic reasons over what’s already known; we make sure what’s known is actually true and traceable to mechanism.”* At the platform layer we are potentially **underneath** them — their agents could in principle query a Graphomics-built federated graph.

-----

## 3. Why the Genentech paper matters (the bridge, not a side credential)

The paper’s three-tier architecture **is** the vertical bridge across the capital/science gap:

- **Tier 1 — Program Decision Layer** = the CEO’s altitude (milestones, Go/No-Go gates).
- **Tier 2 — Assay Protocol Layer** = the CSO’s altitude (mechanism, decision logic, failure modes).
- **Tier 3 — Execution Infrastructure Layer** = the bench reality underneath.

The cross-tier traversal (e.g., `SOURCED_FROM` running from a program milestone down to execution integrity) is a **capital-allocation question answered in root-level scientific terms, automatically.** That is the bridge made operational — deployed in production at Genentech, peer-reviewed. The neat features (e.g. `MASKED_BY` / silent-failure detection) are not the point; they are *proof that the vertical traversal carries real scientific weight* — that rolling science up into a decision doesn’t flatten it into a dashboard metric.

**Honest boundary to hold:** the paper’s bridge was built on deep, expert-elicited assay knowledge at Genentech. The demo’s bridge is built on MicroMap + Reactome + public data — same *architecture*, thinner *substrate*. Claim the architecture generalizes (true). Do **not** imply the public-data demo carries the same evidentiary weight as the Genentech deployment.

Clean way to say both: *“Here’s the architecture proven at depth in production. Here’s the same architecture running live, domain-agnostic, on open science. The engine is identical — point it at your data and it goes as deep as your data allows.”*

-----

## 4. Platform posture (harness-agnostic) + the moat

The platform was engineered so the orchestrator is **swappable** via MCP — Claude Code, Codex, a pharma’s own internal agent, or our own Nexus can all drive the same vertical traversal. This is a *platform* story, not a *product* story, and it’s what earns a platform multiple.

**The moat question a sharp VC will ask:** *“If the orchestrator is swappable, isn’t the orchestrator where the value is — and that’s not you?”*

**The answer (must be ready):** the orchestrator is a commodity; the **mechanistic, federated, provenance-carrying graph is not.** Anyone can point an agent at data. Almost no one can give that agent a substrate where `taxon → metabolite → Reactome pathway → disease` is traversable with edge-level provenance and confidence, federated across domains and institutions. The harness is interchangeable *precisely because* the hard, defensible thing lives in the layer below it (MapForge construction, MicroMap curation, the federation contract). We deliberately commoditized the layer above us.

**Demo consequence:** we drive the live demo with **Claude Code**, not Nexus — to *prove* the platform claim live instead of asserting it. Nexus is mentioned as the purpose-built first-party orchestrator tuned for scientists who don’t live in a terminal. Both stories land in one demo.

-----

## 5. The products in play

- **MapForge** — schema-driven KG engine; builds federated graphs from any source. *This is what we demo building/federating, not MicroMap itself.*
- **MicroMap** — the product *of* MapForge: curated + distributional microbiome KG. The base mechanistic substrate for the demo.
- **Workbench** — computational pipeline execution, invoked via MCP.
- **Nexus** — first-party agentic orchestrator. In scope, secondary for this demo (Claude Code drives).

-----

## 6. Demo design constraints (decided)

1. **No wet lab → literature + omics data are the lab proxy.** Framing: *Graphomics is the in-silico half of a lab-in-the-loop system; the partner’s wet lab plugs into the open end.* Do not imply we run a physical loop.
1. **Domain-agnostic company, single-domain demo.** Show **microbiome** end-to-end (where MicroMap has depth and federation is actually built). Assert generality as an architectural property; do not try to show breadth by showing three shallow domains (that reads as “unfocused” — the exact wedge concern already flagged).
1. **“Mechanistic” must visibly land.** The metabolite-mediated path has to show a **Reactome-grounded** middle hop (`metabolite → named pathway → disease`), not two stacked correlations. If we say mechanistic and show co-occurrence, we lose the room.
1. **Reactome for the demo (CC0, safe).** **KEGG is commercially licensed** — mention as a licensed-tier integration, do **not** embed KEGG-derived data in anything shown to investors until licensing is confirmed.
1. **Provenance heterogeneity is the integrity story, not a liability.** Surface edge-level provenance + confidence type (mined / curated / literature) so the researcher sees which hops are solid pathway biochemistry and which are statistical hypotheses needing validation.

-----

## 7. End-to-end demo flow (microbiome spine)

**Framing question (CEO altitude, crosses the gap):**

> *“We’re considering investing in a microbiome-targeted program for disease X. Is the underlying biology real enough to bet on, and what’s the mechanistic basis?”*

1. **Claude Code (orchestrator)** receives the question — visibly a third-party harness driving our stack.
1. **MapForge / MicroMap (base):** returns ranked taxa + metabolite-mediated path; the `metabolite → pathway → disease` hop **grounded in Reactome**; every edge tagged by provenance type. ← *mechanistic claim lands here.*
1. **Orchestrator (synthesis):** pulls primary literature supporting/contradicting the top hypothesis — literature-as-lab-proxy entering the loop.
1. **Workbench (execution):** runs a real microbiome pipeline (DADA2 → taxonomic profiling → DESeq2 differential abundance) on a real public dataset, testing predicted enrichment in disease vs. control. ← *omics data + literature together are the wet-lab stand-in.*
1. **Orchestrator (interpret):** confirms / refutes / refines against the original mechanistic hypothesis.
1. **MapForge (federation):** the validated finding is written to a **federated source** (the BIOM-derived microbiome network); one live query then spans MicroMap core + the federated instance with unified, two-layer provenance.

**Then roll the answer back up** to the capital-allocation altitude — same traversal serves both the CEO and the CSO. *That* is the trick neither a top-down tool (Perceptic) nor a bottom-up bench tool can do.

-----

## 8. Staging recommendation (live vs pre-staged)

Live KG construction from unstructured data is the slowest, most failure-prone thing we could put on stage. Recommended split:

- **Pre-stage:** federated instance creation (MapForge ingests a source → federated graph exists before the demo).
- **Live:** federation *working* — one query spanning core + federated instance, returning unified results with provenance.
- **Live:** one contribution *write* (step 6) — fast, proves the engine is real, low risk.

**Fallback:** if the full 1→6 loop is flaky a week out, do **not** debut it live. Record the golden path and narrate live with the real product visible. A clean recorded loop beats a broken live one in front of an investor.

-----

## 9. OPEN QUESTIONS — feasibility hunt for cofounders

These determine whether June 9 is a **live** demo or a **recorded golden-path** demo. Ordered by risk. (Claude could not verify these — no repo access in this session.)

### Highest risk — the MCP control surface

- [ ] **Can an *external* harness (Claude Code / Codex) drive the full loop via MCP, or is the control surface Nexus-coupled?** This is the single load-bearing assumption of the whole demo. *(Owner: Varun)*
  - Where to look: the Workbench MCP server definition / tool manifest. Does it expose pipeline invocation as MCP tools any client can call?
- [ ] **Has the full 1→6 loop ever run start-to-finish, unattended, with an external harness driving?** (Known: no, as of now — one week to build. Confirm what specifically is unproven.)

### Workbench / pipelines

- [ ] **Which Workbench pipelines are wired to MCP today?** *(Owner: Varun)*
- [ ] **Do any operate on 16S / shotgun microbiome data (DADA2 / profiling / DESeq2)?** Step 4 is real only if yes. If not, what’s the gap to wire one?

### MicroMap / the “mechanistic” claim

- [ ] **Does the MicroMap query API expose the `taxon → metabolite → disease` traversal?** *(Owner: Demetrius)*
- [ ] **Does each edge carry provenance + confidence type** (mined / curated / literature) retrievably, so we can surface it in the demo?
- [ ] **Is the Reactome `metabolite → pathway → disease` hop ingestible into MapForge for the demo** — already present, or net-new build? *(Owner: Demetrius)*
- [ ] **KEGG licensing:** confirm current commercial terms before any KEGG-derived data goes near the deck. Default to Reactome.

### Federation / the BIOM network

- [ ] **Does the BIOM→network pipeline’s output speak the federation contract** (`:FederatedSource`, capability enum, Neo4j/Bolt shape from the executor contract), or does it need a conversion layer first? *(Owner: colleague building it + Varun)*
  - If it doesn’t speak the contract yet, **this is the integration to prioritize this week** for the step-6 payoff.
- [ ] **What exactly does the orchestrator call to invoke (a) a contribution write and (b) a federated query?** Confirm both are MCP-exposed for an external harness.

### Demo logistics

- [ ] Public microbiome dataset chosen for step 4 (disease X with known disease/control 16S or shotgun data).
- [ ] Decide live-vs-recorded per Section 8 once the above are answered.

-----

## 10. One-line thesis (for the top of the demo + the deck)

> Every pharma has a lossy, expensive translation gap between the people who allocate capital and the people who defend the science. Graphomics is the mechanistic-science substrate that closes it — a federated, tiered, provenance-carrying knowledge layer that lets any agentic harness reason from a CEO-altitude question down to root-level mechanism and back up, without flattening the science. We proved the architecture at depth in production at Genentech. We’re showing it run live, domain-agnostic, on open science, driven by an off-the-shelf agent — because the orchestrator is a commodity and the substrate is the moat.