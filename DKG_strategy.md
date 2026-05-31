# Distributional Knowledge Graph (DKG) Strategy

## Core Thesis

Most analytical systems reduce relationships to a small number of scalar summaries:

* Pearson correlation
* regression slope
* p-value
* mutual information
* differential abundance/expression

This collapses rich biological or system behavior into a weak representation.

The central idea of the Distributional Knowledge Graph (DKG) framework is:

```text id="vhg1md"
relationships themselves are structured objects
```

not scalar coefficients.

A DKG attempts to characterize the full geometry of the relationship between variables:

* mean behavior
* variance structure
* asymmetry
* tail enrichment
* threshold behavior
* regime shifts
* predictive utility
* stability
* distributional divergence
* conditional state occupancy

The result is not merely a correlation network.

It is a graph of typed, distribution-aware relationships.

---

# General Philosophy

## Conventional Pairwise Analysis

Traditional systems ask:

```text id="8g4lh0"
are X and Y associated?
```

The DKG instead asks:

```text id="wb8a4s"
what kind of relationship exists between X and Y?
```

This is a much richer question.

Examples:

* linear monotone effect
* saturation
* threshold activation
* variance expansion
* left-tail vulnerability enrichment
* bifurcation
* multimodality
* subgroup-selective effects
* unstable or context-specific effects

These distinctions are often biologically or mechanistically critical.

---

# The Relationship Representation Concept

For every pair `(X, Y)`:

```text id="vl4o1f"
X -> Y
```

the system computes a relationship phenotype vector.

Example:

```text id="0f2ndz"
[
  pearson_r,
  spearman_rho,
  linear_slope,
  nonlinear_gain,
  variance_ratio,
  tail_risk_difference,
  wasserstein_distance,
  threshold_location,
  predictive_auc,
  bootstrap_stability,
  ...
]
```

This vector becomes the true representation of the edge.

The graph itself is a secondary structure built from these representations.

---

# The Most Important Conceptual Shift

The DKG is not primarily:

```text id="mxfx1k"
a graph of variables
```

It is:

```text id="g5b1zn"
a representation system for relationships
```

This distinction is critical.

The graph is only one projection of the metric space.

The metric space itself is the foundational object.

---

# The 12 Relationship Families

## 1. Marginal Geometry

Characterize individual variables before pairwise analysis.

Purpose:

* QC
* detect pathological distributions
* remove low-information variables
* characterize sparsity and support

Examples:

* entropy
* zero inflation
* skewness
* kurtosis
* dynamic range
* multimodality

---

## 2. Global Linear Association

Classical global relationship summaries.

Examples:

* Pearson
* Spearman
* Kendall
* linear slope
* robust slope

Purpose:

* monotone effect detection
* baseline association structure

---

## 3. Conditional Mean Shape

Characterize nonlinear mean behavior.

Questions:

* saturation?
* sigmoidal response?
* plateau?
* convexity?
* non-monotonicity?

Examples:

* spline improvement
* local slopes
* piecewise fit gain

---

## 4. Conditional Variance Structure

Characterize heteroskedasticity.

Questions:

* does predictor destabilize the system?
* does variance expand in certain regimes?

Examples:

* variance ratio
* variance slope
* residual spread structure

---

## 5. Tail Behavior

Characterize enrichment of extreme states.

Questions:

* does X enrich extreme outcomes?
* does X suppress resistant states?

Examples:

* left-tail enrichment
* percentile-event enrichment
* tail risk ratios
* tail PR-AUC

This is often extremely important in:

* therapeutic vulnerability
* rare events
* subgroup discovery

---

## 6. Skewness and Asymmetry

Characterize directional instability.

Questions:

* do rare severe states emerge?
* is variability directional rather than symmetric?

Examples:

* skewness shift
* quantile asymmetry
* asymmetric spread expansion

---

## 7. Regime and Threshold Structure

Characterize state transitions.

Questions:

* does behavior qualitatively change beyond a threshold?
* is there regime switching?
* does the system saturate?

Examples:

* threshold scans
* piecewise regression
* threshold tail enrichment
* ΔAIC of regime models

---

## 8. Distributional Shift

Compare entire conditional distributions.

Questions:

* does the overall distribution move?
* is there broad redistribution of states?

Examples:

* Wasserstein distance
* KS statistic
* energy distance
* quantile profile shifts

This is often richer than correlation.

---

## 9. Predictive Utility

Measure recoverable information.

Questions:

* can X predict Y?
* does the relationship generalize?

Examples:

* CV-R²
* AUROC
* PR-AUC
* Brier score
* nonlinear prediction gain

---

## 10. Robustness and Stability

Characterize reproducibility of the relationship itself.

Questions:

* does the effect survive perturbation?
* is it driven by outliers?
* is threshold stable?

Examples:

* bootstrap stability
* sign consistency
* threshold variance
* edge recovery frequency

---

## 11. Covariate-Adjusted Relationships

Separate mechanism from confounding.

Questions:

* does the relationship persist after conditioning?
* is the edge merely a lineage proxy?

Examples:

* lineage-adjusted effect retention
* covariate-adjusted predictive gain

---

## 12. Mixture and Cluster Structure

Detect latent state mixtures.

Questions:

* continuous relationship?
* hidden subgroups?
* bimodality?

Examples:

* multimodality metrics
* cluster-separation geometry
* within-state vs between-state effects

---

# Graph Construction Philosophy

The DKG should NOT initially construct a graph.

The correct order is:

```text id="7wm0kk"
1. characterize relationships
2. store relationship vectors
3. analyze metric structure
4. define edge ontology
5. admit graph edges
6. construct graph projections
```

This is fundamentally different from:

* correlation networks
* coexpression graphs
* nearest-neighbor graphs

---

# Typed Relationship Ontology

The DKG enables explicit edge classes.

Examples:

* linear_edge
* threshold_edge
* saturating_edge
* variance_expanding_edge
* tail_enrichment_edge
* bifurcation_edge
* predictive_edge
* unstable_edge
* lineage_proxy_edge

This ontology becomes extremely important for:

* querying
* prioritization
* graph algorithms
* downstream ML
* mechanistic interpretation

---

# Relationship Embedding Space

A major extension of the DKG concept is:

```text id="mjlwmx"
relationships become learnable objects
```

The edge vectors themselves can be:

* clustered
* embedded
* classified
* predicted

Examples:

* UMAP over edge vectors
* edge archetype discovery
* edge retrieval systems
* edge embeddings

This transforms the DKG from:

```text id="4jlwmx"
a statistics engine
```

into:

```text id="5jlwmx"
a relationship representation system
```

---

# Recursive Relationship Learning

One of the deepest extensions of the DKG idea is:

```text id="6jlwmx"
relationships can predict relationships
```

Examples:

* expression-expression structure predicts expression-response structure
* edge neighborhoods predict missing edge types
* graph geometry predicts mechanistic classes

This moves toward:

* self-supervised relationship learning
* transferable relationship embeddings
* graph representation learning

---

# Adaptive Compute Allocation

The DKG should not compute expensive analyses uniformly.

Instead:

```text id="7jlwmx"
cheap metrics nominate expensive analyses
```

Example:

* Pearson screen
* then threshold scan
* then GAMLSS refinement
* then interaction analysis

This creates:

* scalable computation
* resource-aware analysis
* agentic analytical workflows

---

# Domain Generality

The DKG is intentionally domain-agnostic.

Applicable domains include:

* cancer dependency biology
* microbiome studies
* metabolomics
* ecology
* finance
* recommender systems
* game telemetry
* MTG draft systems
* clinical response modeling

The framework is fundamentally about:

```text id="8jlwmx"
distribution-aware relational structure
```

not any specific biological modality.

---

# DepMap-Specific Strategy

## Core Objective

Construct a mechanism-aware vulnerability atlas linking:

```text id="9jlwmx"
molecular state
→ selective dependency
```

rather than merely:

```text id="a0jlwmx"
expression
→ correlation
```

---

# Most Valuable Relationship Type

The highest-value edges are:

```text id="b1jlwmx"
expression -> Chronos
```

These directly support:

* target discovery
* biomarker discovery
* patient stratification
* synthetic lethal inference

---

# Expression -> Expression Layer

This layer is primarily:

* mechanistic
* contextual
* aggregative

Purpose:

* define biological programs
* support mechanism hypotheses
* denoise single-gene signals
* identify compensatory structure

---

# Therapeutically Important Edge Archetypes

## 1. Subgroup-Selective Vulnerability

Example:

```text id="c2jlwmx"
high PPARG
→ RXRA/RXRB vulnerability
```

Characteristics:

* strong tail enrichment
* subgroup separation
* moderate global slope
* stable distributional shift

Often more actionable than:

```text id="d3jlwmx"
globally strong diffuse effects
```

---

## 2. Paralog Buffering

Example:

```text id="e4jlwmx"
low paralog A
→ dependency on paralog B
```

Key signals:

* one-sided enrichment
* threshold behavior
* percentile-event predictability
* stable left-tail structure

---

## 3. State-Defined Vulnerability Programs

Use expression-expression edges to define:

* EMT states
* stress programs
* lineage programs
* metabolic programs

Then attach:

```text id="f5jlwmx"
program -> dependency
```

edges.

This is often more robust than single genes.

---

# Key Prioritization Principle

Do NOT prioritize:

* strongest Pearson
* strongest mean shift
* lowest p-value

Prioritize:

```text id="g6jlwmx"
stable, subgroup-selective,
mechanistically coherent vulnerability geometry
```

---

# Recommended DepMap Tables

## 1. `relationship_vectors`

One row per edge.

Contains:

* all relationship metrics
* robustness summaries
* edge ontology labels

---

## 2. `edge_embeddings`

Learned or compressed edge representations.

Used for:

* clustering
* retrieval
* edge archetype analysis

---

## 3. `mechanism_context`

Mechanistic support layer:

* nearby programs
* pathway structure
* module membership
* lineage context

---

## 4. `therapeutic_hypotheses`

One row per:

```text id="h7jlwmx"
biomarker state -> target dependency
```

Contains:

* responder subgroup
* dependency geometry
* mechanism support
* predictive utility
* robustness
* tractability

---

# Long-Term Vision

The long-term vision is NOT merely:

```text id="i8jlwmx"
a better correlation graph
```

It is:

```text id="j9jlwmx"
a universal representation framework
for structured relationships
```

where:

* relationships are typed objects
* graphs emerge from relationship geometry
* edge embeddings become transferable
* mechanistic reasoning operates over relationship structure itself

The therapeutic vulnerability atlas is one application of this broader framework.
