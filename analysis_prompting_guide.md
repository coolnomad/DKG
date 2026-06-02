# AI-Assisted Biological Interpretation Guide

Companion to `analysis_protocol.md`. Describes how to prompt an AI agent to extract
mechanistic and biological insight from DKG outputs at each analysis step.

The NMT1 analysis is the reference case throughout. Prompts are written generically
but include NMT1 examples to show what good output looks like.

---

## How to use this guide

Run the mechanical steps in `analysis_protocol.md` first. Then use the prompts here
to interpret what you got. The most productive pattern is:

1. Show the agent the raw numbers (top hits, key metrics)
2. Ask a specific mechanistic question — not "what does this mean" but "why would these
   two genes be co-essential" or "what's the mechanism by which X could sensitize to Y"
3. Follow up on anything surprising or that contradicts your prior
4. Ask the agent to log findings before moving to the next step

The agent has no memory of prior sessions unless you paste context or point it at the
session log. Always orient it at the start of a new session.

---

## Session startup prompt

Paste this (filled in) at the start of any new session:

```
We are doing a DKG (distributional knowledge graph) analysis of [TARGET] dependency
in DepMap. The pipeline screens expression features (XY mode) and dependency features
(YY' mode) against [TARGET] CRISPR chronos scores across cell lines.

Key files:
- output/[TARGET]_full/tier2_target_full.parquet — XY Tier 2 shape metrics
- output/[TARGET]_yy_full/tier2_target_full.parquet — YY' Tier 2 shape metrics
- session_log_[TARGET].md — running log of findings

Prior findings summary: [paste 3-5 bullet points from the log]

Today I want to [goal].
```

---

## Step 1 — Interpreting XY Tier 2 results

### Orienting prompts

After reading the top positive and negative Pearson correlates:

> "The top positive correlates of [TARGET] dependency are [genes]. The top negative
> correlates are [genes]. What biological story do these tell about when [TARGET] is
> essential vs dispensable?"

> "We see [GENE] as the top positive correlate (r=+X, wasserstein=+Y, delta_aic=Z).
> Walk me through what each of those three metrics is telling us about the nature of
> this relationship."

> "[GENE] has r=+0.4 but delta_aic=15. What does a high delta_aic mean for patient
> selection strategy compared to a gene with r=+0.4 but delta_aic=0?"

### Biomarker type interpretation

> "The wasserstein shift for [GENE] is negative (wass=-0.6). What does that mean in
> terms of how I'd use this as a patient selection biomarker?"

> "We have two classes of correlates: one cluster with high wasserstein shift and low
> IQR ratio, another with low wasserstein and high IQR ratio. What do these two classes
> represent mechanistically, and which is more actionable for patient selection?"

### Surprises and contradictions

> "[GENE] is a positive correlate of [TARGET] dependency but I would have expected it
> to be a negative correlate based on [prior]. Why might the data show this?"

> "The top correlates don't include [EXPECTED_GENE] which is the canonical biomarker
> for this context. Where does [EXPECTED_GENE] rank and what might explain its absence
> from the top hits?"

---

## Step 2 — Interpreting XY community structure

After running Louvain community detection on XY biomarker genes:

> "We have [N] communities from the XY analysis. Community [C] contains [genes] and
> Community [D] contains [genes]. Are these biologically distinct programs or
> different readouts of the same thing?"

> "Community [C] enriches for [pathway from Enrichr]. Community [D] enriches for
> [pathway]. What does it mean for [TARGET] biology that both programs independently
> predict its essentiality?"

> "Which of these communities is most actionable for patient selection and why?"

---

## Step 3 — Interpreting XX shape communities

After running `xx_shape_communities.py --mode genes` and `--mode pairs`:

> "The gene-level shape analysis puts [GENE_A] and [GENE_B] in the same community
> (wass=X, iqr_ratio=Y). But the pair-level analysis shows [GENE_A] x [GENE_B] in a
> different cluster than [GENE_A] x [GENE_C]. What does this hub heterogeneity tell
> us about [GENE_A]'s biology?"

> "In the pair-level shape communities, cluster [N] has high delta_aic (=16) and
> contains [ZEB1 x LIX1L, ZEB1 x IKBIP]. What does it mean that ZEB1's relationship
> with these genes is threshold-driven rather than linear?"

> "We found [TARGET]'s co-expression partners fall into [N] shape communities.
> Which shape pattern is most consistent with a mechanism where [TARGET] loss is
> catastrophic only after EMT commitment / only in high-proliferation cells / etc.?"

---

## Step 4 — Interpreting YY' co-essentiality results

After reading top positive and negative correlates from `output/TARGET_yy_full/`:

### Co-essential genes (positive correlates)

> "The top co-essential genes with [TARGET] are [list]. What pathway or complex do
> these point to? Is there a coherent biological story?"

> "[GENE] is co-essential with [TARGET] at r=+X. Are these paralogs? If so, what
> does their co-essentiality tell us — are they in the same complex or do they
> provide redundant functions in different contexts?"

> "We see [KLF_FAMILY] genes dominating the co-essential list. What does KLF family
> co-essentiality tell us about the cell identity context where [TARGET] is essential?"

### Anti-correlated genes (sensitizer candidates)

> "The top anti-correlated dependency genes are [list]. What biological program do
> these represent, and why would cells dependent on this program NOT need [TARGET]?"

> "[GENE_A] anti-correlates with [TARGET] at r=-X with delta_aic=Y. The high delta_aic
> suggests a threshold. What does that threshold structure mean for how we'd design
> a combination treatment using a [GENE_A] inhibitor?"

> "We have [GENE_A] (r=-X, d_aic=15) and [GENE_B] (r=-X, d_aic=1) both as
> anti-correlates. Why does the threshold structure make [GENE_A] a better combination
> candidate than [GENE_B] even if their correlations are similar?"

### Combination strategy

> "[GENE_A] and [GENE_B] both anti-correlate with [TARGET]. Are they the same pathway
> at different levels, or independent escape mechanisms? Does it make more sense to
> combine [TARGET] inhibitor with [GENE_A] inhibitor, [GENE_B] inhibitor, or both?"

> "Explain the mechanism by which a [GENE_A] inhibitor could sensitize cells that are
> currently [TARGET]-independent to become [TARGET]-dependent. What is the cell
> biological chain of events?"

> "What existing clinical-stage inhibitors exist for [top anti-correlate genes]?
> Which has the best safety profile for a combination trial?"

---

## Step 5 — Interpreting TCGA prevalence results

After running `tcga_marker_prevalence.py`:

> "The top TCGA cohorts by [MARKER]-low prevalence are [list with %]. Which of these
> are mechanistically plausible vs likely to be false positives based on the marker
> biology? Flag any that are likely lineage biology rather than the mechanism we care about."

> "[COHORT] has [X]% [MARKER]-low but VIM and ZEB1 are both elevated there. Does
> that support or contradict the [mechanism] model?"

> "We want to prioritize indications by NMT2-low prevalence AND mechanistic plausibility
> AND unmet medical need AND cohort size. Given [data], rank the top 5 indications."

### Subtype stratification

> "The bulk [COHORT] number is [X]%. This masks subtype heterogeneity. Which subtypes
> within [COHORT] are most likely to have elevated [MARKER]-low prevalence based on
> the biology, and how would I go about checking this in cBioPortal?"

> "We found [SUBTYPE_A] has [X]% [MARKER]-low but [SUBTYPE_B] only has [Y]%. This is
> the opposite of what I expected. What biological explanation reconciles this with
> what we know about [SUBTYPE_A] biology?"

---

## Step 6 — Patient selection and indication strategy

After all analyses are complete:

> "Synthesize everything we've found. What is the primary patient selection biomarker
> for [TARGET] inhibition? What are the top 3 indications? What is the combination
> strategy for patients who don't fit the primary biomarker?"

> "We have two non-overlapping patient populations: [BIOMARKER_A]-selected (mechanism
> 1) and [BIOMARKER_B]-selected (mechanism 2). Are these actually complementary
> indications, or could they be run as a single biomarker-stratified trial?"

> "What would a phase 1/2 trial design look like for [TARGET] inhibitor in [INDICATION]
> using [BIOMARKER] as the selection criterion? What is the key pharmacodynamic readout?"

> "Log this patient selection strategy to the session log."

---

## Additional queries recommended based on NMT1 learnings

These are analyses we did not run for NMT1 but would add value for any target:

### Druggability of co-essential and anti-correlated genes

> "Of the top 20 co-essential genes with [TARGET], which are druggable — either by
> existing approved drugs, clinical-stage compounds, or tool compounds? Rank by
> combination attractiveness."

> "Of the top anti-correlates, which have FDA-approved inhibitors already in oncology
> trials? Focus on ones with delta_aic > 5."

### Paralog analysis

> "What are the known paralogs of [TARGET]? Are any of them in the co-essential or
> anti-correlated YY' lists? What does their position tell us about the synthetic
> lethality mechanism?"

> "Is [TARGET] itself in any paralog pair columns in the DepMap CRISPR screen
> (e.g. TARGET_PARALOG)? If so, what is the co-essentiality of the pair score vs
> the individual gene score?"

### Mechanistic hypotheses for experimental validation

> "Based on everything in this analysis, what are the top 3 mechanistic hypotheses
> that explain why [TARGET] is essential in [BIOMARKER_STATE] cells? State each as
> a testable hypothesis with a specific experiment."

> "What loss-of-function or gain-of-function experiment in a cell line model would
> most efficiently validate the [BIOMARKER] → [TARGET] essentiality relationship?"

### Cell line model selection

> "Given the biomarker profile we've identified ([BIOMARKER_A]-low, [MARKER_B]-high),
> which DepMap cell lines best represent the target patient population for in vitro
> validation? What are the top 5 and what are their caveats?"

> "Are there publicly available PDX or organoid models that match the [INDICATION] +
> [BIOMARKER] profile we've identified?"

### Cross-target comparison

> "We previously ran this analysis for [OTHER_TARGET]. How do the XY biomarker
> communities compare between [TARGET] and [OTHER_TARGET]? Do they share co-essential
> genes or sensitizer candidates?"

> "Could [TARGET] and [OTHER_TARGET] be co-targeted in the same indication? What
> would need to be true about their biomarker populations for this to work?"

### Shape metric deep dives

> "Pull the tier2 row for the [GENE_A] x [TARGET] pair and walk me through every
> shape phase metric that has a large value. Build a mechanistic narrative for what
> this specific relationship looks like across cell lines."

> "Which pairs in our tier2 output have the highest delta_aic? For each, interpret
> what the threshold means — is it a cell state transition, an expression threshold,
> or a co-dependency threshold?"

### Session log hygiene

> "Read the session log and tell me if any earlier conclusions should be revised or
> flagged as superseded based on what we've found since."

> "Summarize the session log into a 10-bullet executive summary suitable for sharing
> with a collaborator who hasn't seen this analysis."

---

## Logging discipline

Ask the agent to log after every major finding, not at the end of the session.
Good prompts:

- "Log this finding to the session log."
- "Log this and update the indication priority if it changes."
- "Log the mechanistic interpretation we just discussed."

The session log is the primary artifact. Code outputs are reproducible; the biological
interpretation is not. If it's not in the log, it's lost.

---

## What makes a good vs. poor prompting session

| Good | Poor |
|------|------|
| Show raw numbers, ask specific question | Ask "what does this mean" without data |
| Follow up on surprises and contradictions | Accept first answer without probing |
| Ask for mechanism, not just description | Ask agent to "summarize" without direction |
| Ask agent to flag false positives | Assume all statistical hits are real |
| Log before moving on | Reconstruct at the end from memory |
| Ask what experiments would test the hypothesis | Stop at the hypothesis |
| Ask about the prior literature for key hits | Treat computational findings as standalone |
