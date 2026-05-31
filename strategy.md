# Strategy for Using Pairwise Expression and Chronos Metrics

## Goal

Use the stored pairwise metrics from:

- expression -> Chronos
- expression -> expression

to build a mechanism-aware vulnerability atlas that can identify:

- potential drug targets
- biomarker-defined responder subgroups
- candidate synthetic lethal relationships
- resistance mechanisms
- mechanistic explanations for observed dependencies

The key idea is that this framework is more powerful than a correlation screen because it captures:

- linear effect size
- nonlinearity
- threshold behavior
- tail enrichment
- distributional shift
- predictive utility
- stability under resampling

This makes it possible to prioritize subgroup-selective therapeutic vulnerabilities rather than just globally associated pairs.

## Core Principle

The most important use of this analysis is not to ask:

"Which genes are correlated?"

It is to ask:

"Which molecular states define a subset of models that become selectively dependent on a targetable gene?"

That is the right abstraction for drug target discovery in DepMap-like systems.

## Best Use Case

The highest-value product is a biomarker-to-vulnerability map built from expression -> Chronos edges.

For each predictor-target pair, the practical question is:

- does an expression state predict selective essentiality of a gene?
- is the effect thresholded, nonlinear, or one-sided?
- does it identify a clinically useful subgroup rather than a weak global average trend?
- is the signal stable and predictive out of sample?
- is there mechanistic support from the expression-expression network?

A strong hit is not merely significant. A strong hit has:

- a meaningful dependency shift
- a clear biomarker regime
- robust predictive performance
- resampling stability
- a plausible mechanism

## Most Powerful Discovery Modes

### 1. Biomarker-Defined Target Discovery

This is the most direct route to target plus patient-selection biomarker.

Look for targets whose Chronos score becomes strongly negative only in a defined expression state.

Useful metrics:

- `left_tail_risk_difference`
- `left_tail_risk_ratio`
- `wasserstein_1`
- `mean_shift`
- `regime_threshold`
- `regime_delta_aic`
- `left_tail_auc`
- `left_tail_pr_auc`
- stability summaries from phase 10

Interpretation:

- A large mean shift alone is not enough.
- A thresholded subgroup with strong left-tail enrichment is often more therapeutically meaningful.
- A pair with moderate correlation but sharp event enrichment can be more actionable than a pair with stronger global correlation.

### 2. Synthetic Lethal and Paralog Buffering Discovery

This is one of the most compelling applications.

Examples:

- low expression of gene A induces dependency on gene B
- loss-like expression state of one paralog creates dependency on the other
- compensatory pathway activation creates a targetable liability

Look for:

- nonlinear or thresholded activation of dependency
- strong left-tail enrichment in the vulnerable state
- one-sided asymmetry rather than symmetric shifts
- stability across bootstrap and subsampling

Expression-expression edges can then be used to test whether the surrounding transcriptional context supports a compensation mechanism.

### 3. State-Specific Vulnerability Programs

Single-gene biomarkers are useful, but many real biological dependencies are state-based.

Expression-expression edges can define:

- transcriptional programs
- pathway modules
- stress states
- lineage-like states

Then expression -> Chronos edges can be used to attach selective vulnerabilities to those states.

This makes it possible to discover patterns like:

- EMT-like state -> dependency on receptor or kinase X
- interferon-high state -> dependency on immune-stress pathway Y
- oxidative phosphorylation program -> dependency on mitochondrial target Z

This is often more biologically meaningful and more translatable than single-gene associations.

### 4. Resistance and Escape Mechanism Inference

The same atlas can be used to understand why a dependency disappears in a subset of models.

For a known target, ask:

- which expression states weaken dependency?
- which alternate pathways are active in those states?
- is there a regime boundary that separates sensitive and resistant models?

This is useful for:

- resistance mechanism discovery
- combination target nomination
- biomarker-aware patient segmentation

### 5. Therapeutic Geometry-Based Target Prioritization

Not every statistically detectable association is useful.

The best therapeutic targets often show:

- a clean responder subgroup
- a biomarker that separates responders from non-responders
- strong selective dependency within that subgroup
- relative lack of broad dependency outside that subgroup
- mechanistic coherence with the surrounding network

This framework is useful because it can distinguish:

- diffuse weak trends
- nonlinear but real subgroup effects
- thresholded therapeutic regimes

## Recommended Data Product

The most powerful end product is a heterogeneous graph or edge table system with typed relationships.

### Nodes

- genes
- dependency targets
- expression modules or pathway programs
- lineages
- later: mutations, copy number states, fusions, drug response features

### Edge Types

- expression -> Chronos
- expression -> expression

### Edge Attributes

Store the important metrics for each edge:

- slope and sign
- rank-based association
- nonlinear fit improvement
- threshold location
- tail risk metrics
- distribution shift metrics
- predictive performance metrics
- stability metrics

This gives a mechanism-aware vulnerability graph rather than an unstructured results dump.

## Practical Analysis Strategy

### Step 1. Run the Pairwise Screen at Scale

Materialize pairwise metrics for:

- expression -> Chronos across all genes or a focused target set
- expression -> expression across all genes or a biologically constrained set

At minimum, store:

- phase 2 association metrics
- phase 3 nonlinear mean-shape metrics
- phase 5 tail metrics
- phase 7 regime metrics
- phase 8 distribution-shift metrics
- phase 9 predictive metrics
- phase 10 stability summaries

### Step 2. Rank Expression -> Chronos Edges for Target Potential

Do not rank by Pearson correlation alone.

Build a composite prioritization score emphasizing:

- dependency selectivity
- subgroup enrichment
- predictive power
- threshold sharpness
- robustness

An effective target-priority score should reward:

- large negative dependency shifts in one molecular regime
- high left-tail event enrichment
- strong out-of-sample event prediction
- stable effect direction and magnitude
- evidence that the dependency is not just a diffuse lineage effect

### Step 3. Use Expression -> Expression Edges for Mechanism

Once a biomarker-target edge is interesting, use the expression-expression network to ask:

- what genes track with the biomarker?
- what pathways or programs define the same state?
- does the local network support compensation, lineage, stress response, or pathway rewiring?

This converts a vulnerability edge into a mechanism hypothesis.

The main value of the expression-expression layer is explanatory and aggregative:

- explain why the target dependency appears
- define broader state modules
- denoise single-gene signals into pathway-level biology

### Step 4. Collapse Single Genes into Modules

Single-gene expression is noisy and can be biologically narrow.

A strong extension is:

1. use expression-expression edges to infer modules or programs
2. compute module scores per model
3. rerun module -> Chronos analyses

This will often produce:

- stronger predictive signals
- cleaner mechanistic interpretation
- better portability across cohorts

### Step 5. Separate Mechanism from Lineage

This is critical.

Many top expression-dependency hits in DepMap will be lineage-driven rather than mechanistically informative.

At a minimum, compare:

- raw analyses
- lineage-adjusted analyses

Recommended adjustments may include:

- lineage
- culture type
- copy number burden
- proliferation-related covariates if relevant

The lineage-adjusted results are often the most useful for discovering generalizable target-mechanism relationships.

### Step 6. Focus on One-Sided Vulnerabilities

Many clinically useful dependencies are one-sided:

- the vulnerable state is highly dependent
- the opposite state is not especially informative

This makes the tail and regime metrics especially valuable.

In practice, prioritize cases where:

- a defined biomarker state sharply enriches extreme dependency
- the effect is not merely a mild symmetric trend
- the subgroup is large enough to matter biologically and therapeutically

## What a Strong Hit Looks Like

A compelling hit often has the following structure:

- biomarker state: low or high expression of gene A, or a module score regime
- target: gene B Chronos dependency
- evidence:
  - moderate or strong association
  - meaningful nonlinear or thresholded structure
  - strong left-tail enrichment in the vulnerable subgroup
  - good cross-validated predictive performance
  - stable bootstrap and subsampling behavior
  - mechanistic support from nearby expression-expression edges

Example archetypes:

- low `SMARCA4` expression -> `SMARCA2` dependency
- low paralog A expression -> paralog B dependency
- mesenchymal program high -> selective kinase dependency

Each of these suggests:

- a drug target
- a responder biomarker
- a mechanism hypothesis
- a path toward validation

## What Not to Do

- Do not prioritize by Pearson correlation alone.
- Do not confuse lineage markers with causal mechanisms.
- Do not trust unstable effects.
- Do not stop at pairwise ranking without mechanism annotation.
- Do not assume the best target is the strongest average dependency; subgroup-selective effects are often more valuable.

## Suggested Output Tables

To make this system useful at scale, materialize a small number of queryable tables.

### 1. `expr_to_chronos_edges`

One row per expression predictor and Chronos target pair.

Recommended fields:

- predictor gene
- target gene
- association metrics
- nonlinear metrics
- tail metrics
- regime metrics
- distribution-shift metrics
- predictive metrics
- stability summaries
- raw and lineage-adjusted versions

### 2. `expr_to_expr_edges`

One row per expression-expression pair.

Recommended fields:

- predictor gene
- target gene
- association metrics
- nonlinear metrics
- regime metrics
- stability summaries

### 3. `target_priority_table`

One row per candidate therapeutic hypothesis.

Recommended fields:

- target gene
- biomarker gene or module
- direction of biomarker state
- target-priority composite score
- subgroup size
- dependency effect size
- event enrichment
- predictive performance
- stability score
- lineage-adjusted status
- target tractability annotation

### 4. `mechanism_annotation_table`

One row per nominated target-biomarker relationship with mechanism context.

Recommended fields:

- biomarker
- target
- nearby expression program
- candidate pathway
- synthetic lethal / buffering / lineage / resistance interpretation
- supporting expression-expression neighbors

## Ranking Framework for Drug-Target Potential

For target discovery, the ranking objective should favor subgroup-selective vulnerabilities with mechanistic support.

A practical ranking score can combine:

- dependency severity in the predicted responder state
- tail-event enrichment
- distributional separation
- predictive AUC or PR-AUC
- threshold clarity
- resampling stability
- lineage-adjusted persistence
- tractability of the target

Conceptually:

`target_score = selective_dependency + predictiveness + stability + mechanism_support + tractability - confounding_penalty`

This does not need to be perfect initially. The main requirement is that it not collapse back to simple correlation ranking.

## End Product

The strongest end product is a ranked catalog of biomarker-defined therapeutic opportunities, where each entry contains:

- target gene
- biomarker state
- expected responder subgroup
- dependency geometry
- predictive performance
- mechanism annotation
- confidence and stability
- possible resistance state
- possible combination partner
- tractability or existing druggability context

That turns the pairwise screen into a target discovery engine rather than a descriptive analysis.

## Recommended Next Steps

1. Standardize the on-disk schema for pairwise results.
2. Run a pilot on a biologically focused subset before scaling to all pairs.
3. Add lineage-adjusted analyses early.
4. Build a first-pass target-priority score for expression -> Chronos edges.
5. Add mechanism annotation from expression-expression neighborhoods.
6. Collapse single genes into modules for a second-pass analysis.
7. Validate top hits against known synthetic lethal and biomarker-defined target examples.
