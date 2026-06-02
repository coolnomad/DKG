# DKG Target Analysis Protocol

Replication guide for the full single-target distributional knowledge graph analysis.
Developed using NMT1 as the reference case. Replace `TARGET` / `ENTREZ` throughout.

---

## Prerequisites

```bash
# Activate venv
.venv/Scripts/activate        # Windows
source .venv/bin/activate     # Unix

# Key data files
data/processed/xp_filtered.feather        # expression matrix (277 x 20,303)
data/processed/chronos_filtered.feather   # dependency matrix, expression-intersected (277 x 11,745)
data/processed/chronos_26Q1_full.feather  # full dependency matrix (1,538 x 18,532)
```

To (re)build `chronos_26Q1_full.feather` from the raw DepMap release:

```python
import polars as pl
df = pl.read_ipc('C:/GitHub/DepMap/data/26Q1/CRISPR_26Q1.feather')
new_cols = [c.split('..')[0] if '..' in c else c for c in df.columns]
df = df.rename(dict(zip(df.columns, new_cols))).rename({'row_id': 'ModelID'})
df.write_ipc('data/processed/chronos_26Q1_full.feather')
```

---

## Step 1 — XY Target Mode (expression → dependency)

Screens which expression features predict the target's CRISPR dependency score.
Produces per-pair distributional shape metrics (Tier 2) across 5 CV folds + full data.

```bash
python -m dkg \
    --mode target \
    --x-matrix data/processed/xp_filtered.feather \
    --y-matrix data/processed/chronos_filtered.feather \
    --target-col TARGET \
    --output-dir output/TARGET_full \
    --skip-tier3
```

**Key outputs:**
- `output/TARGET_full/tier2_target_full.parquet` — all predictors, full-data Tier 2 shape metrics
- `output/TARGET_full/tier2_target_fold{0-4}.parquet` — CV fold results
- `output/TARGET_full/splits.parquet` — fold assignments

**Reading results (pearson_r is nested in p2 struct):**

```python
import polars as pl
df = pl.read_parquet('output/TARGET_full/tier2_target_full.parquet')
p2 = df['p2_symmetric_pair_metrics'].struct.unnest()
df = df.drop('p2_symmetric_pair_metrics').hstack(p2)
# Sort by |r|, inspect p8_signed_wasserstein_shift, p7_delta_aic
```

**Key Tier 2 fields to inspect:**
- `pearson_r` — linear correlation (positive = presence biomarker, negative = absence biomarker)
- `p8_signed_wasserstein_shift` — distributional shift direction/magnitude
- `p7_delta_aic` — evidence for threshold/piecewise structure (>5 = meaningful)
- `p4_iqr_ratio_high_low` — variance modulation (>1.5 = high-X expands variance)

---

## Step 2 — XY Community Detection (Louvain on biomarker genes)

Groups the top XY biomarker genes into communities by co-expression structure.
Script: `scripts/nmt1_biomarker_communities.py` (NMT1-specific; adapt for new targets).

Alternatively, use the enrichment approach on tier2 output directly (Step 3).

**Inspect community membership:**

```python
df = pl.read_parquet('output/TARGET_full/biomarker_communities.parquet')
# Columns: gene, community_id, pearson_r, p8_signed_wasserstein_shift, ...
```

---

## Step 3 — Enrichr Pathway Enrichment on XY Communities

Script: `scripts/nmt1_enrichr.py` (adapt gene list extraction for new targets).

Submits gene lists per community to Enrichr API, returns enriched gene sets across
MSigDB Hallmarks, KEGG, GO Biological Process, Reactome.

```bash
python scripts/nmt1_enrichr.py \
    --tier2 output/TARGET_full/tier2_target_full.parquet \
    --output-dir output/TARGET_full/enrichr
```

**Key libraries to check:** HALLMARK, KEGG_2021_Human, GO_Biological_Process_2023,
Reactome_2022. Look for coherent pathway enrichment per community.

---

## Step 4 — XX Mode (within-matrix co-expression on biomarker gene subsets)

Run Tier 2 co-expression analysis on the XY biomarker community gene subsets.
Reveals which co-expression structure underlies each XY community.

First, extract the gene list for each community and create a subset feather:

```python
import polars as pl
# Filter expression matrix to genes in community C (e.g. C3+C4)
genes = [...]  # from biomarker_communities.parquet
xp = pl.read_ipc('data/processed/xp_filtered.feather')
subset = xp.select(['ModelID'] + [g for g in genes if g in xp.columns])
subset.write_ipc('output/TARGET_full/xx_c34_input.feather')
```

```bash
python -m dkg \
    --mode xx \
    --x-matrix output/TARGET_full/xx_c34_input.feather \
    --output-dir output/TARGET_full/xx_c34 \
    --skip-tier3
```

**Key outputs:**
- `output/TARGET_full/xx_c34/tier2_deep.parquet` — all pairs with full shape metrics
- `output/TARGET_full/xx_c34/shape_communities.parquet` — gene-level shape communities

**Flags:**
- `--skip-tier3` — skip bootstrap (recommended for exploration)
- `--tier1-threshold 0.0` — include all pairs regardless of |r| (for shape analysis on small gene sets)

---

## Step 5 — XX Shape-Based Community Detection

Clusters genes (or pairs) by their distributional shape profile across 66 Tier 2
shape features. Two modes: aggregate per gene (fast) or per pair (more granular).

```bash
# Gene-level (aggregate shape per gene, then cluster)
python scripts/xx_shape_communities.py \
    --tier2 output/TARGET_full/xx_c34/tier2_deep.parquet \
    --output-dir output/TARGET_full/xx_c34

# Pair-level (each pair as a node — reveals hub heterogeneity)
python scripts/xx_shape_communities.py \
    --tier2 output/TARGET_full/xx_c34/tier2_deep.parquet \
    --output-dir output/TARGET_full/xx_c34 \
    --mode pairs

# Combined multiple subsets
python scripts/xx_shape_communities.py \
    --tier2 output/TARGET_full/xx_c1/tier2_deep.parquet \
            output/TARGET_full/xx_c34/tier2_deep.parquet \
    --output-dir output/TARGET_full/xx_combined_shape
```

**Key parameters:**
- `--sim-threshold 0.5` — cosine similarity cutoff for graph edges
- `--resolution 1.0` — Louvain resolution (higher = more communities)
- `--mode genes` (default) vs `--mode pairs`

**Interpreting communities by shape metrics:**
- High `wass` (p8_signed_wasserstein_shift) → distributional shift-dominated
- High `iqr_ratio` (p4_iqr_ratio_high_low) → variance-modulating relationships
- High `delta_aic` (p7_delta_aic) → threshold/piecewise structure
- Low all metrics → linear co-expression, no distributional complexity

---

## Step 6 — Joint XX+XY Graph

Builds a single NetworkX graph combining XX co-expression edges (positive, |r|≥0.6)
and XY dependency edges (|r|≥0.2), then runs Louvain community detection.

```bash
python scripts/nmt1_joint_graph.py \
    --xx-parquet  output/TARGET_full/tier1_expression.parquet \
    --xy-parquet  output/TARGET_full/tier2_target_full.parquet \
    --target-col  "TARGET..ENTREZ." \
    --output-dir  output/TARGET_full/joint_graph \
    --edge-threshold 0.6 \
    --xy-threshold   0.2
```

**Key outputs:**
- `output/TARGET_full/joint_graph/communities.parquet` — node community assignments
- `output/TARGET_full/joint_graph/joint_edges.parquet` — all edges with type (xx/xy)

**Reading the graph:**
```python
df = pl.read_parquet('output/TARGET_full/joint_graph/communities.parquet')
# target_node = row where gene == TARGET
# Check: degree, betweenness_centrality, community_id
# Check XY neighbors: which XX communities do they land in?
```

---

## Step 7 — YY' Mode (dependency → dependency)

Screens all CRISPR dependency genes as predictors of the target's dependency score.
Reveals co-essential genes (positive r) and sensitizer candidates (negative r).

```bash
python -m dkg \
    --mode target \
    --x-matrix data/processed/chronos_26Q1_full.feather \
    --y-matrix data/processed/chronos_26Q1_full.feather \
    --target-col TARGET \
    --output-dir output/TARGET_yy_full \
    --skip-tier3
```

**Reading results — filter to high-quality pairs:**

```python
df = pl.read_parquet('output/TARGET_yy_full/tier2_target_full.parquet')
p2 = df['p2_symmetric_pair_metrics'].struct.unnest()
df = df.drop('p2_symmetric_pair_metrics').hstack(p2)

good = df.filter(
    (pl.col('x_col') != 'TARGET') &
    (pl.col('p3_n') >= 200) &
    (pl.col('p8_signed_wasserstein_shift').is_not_null())
)
# Sort by pearson_r descending (co-essential) or ascending (anti-correlated/sensitizers)
```

**Interpreting results:**
- **Positive correlates** — co-essential genes; same pathway/complex as target
- **Negative correlates** — anti-correlated; cells dependent on these are NOT dependent
  on target. Druggable anti-correlates = sensitizer combination candidates.
- **High d_aic on anti-correlates** — threshold structure; sensitizer doesn't need full
  target engagement, just needs to push cells past a regime boundary.

**Note:** Paralog pair columns (e.g. `TARGET_PARALOG`) represent combined CRISPR KO
scores in DepMap, not individual gene scores.

---

## Step 8 — TCGA Pan-Cancer Marker Prevalence

Fetches pan-cancer z-score expression from cBioPortal for any gene set and computes
% low / % high / median z per TCGA PanCan Atlas cohort (30 cohorts, ~10K samples).

```bash
python scripts/tcga_marker_prevalence.py \
    --genes TARGET:ENTREZ GENE2:ENTREZ2 GENE3:ENTREZ3 \
    --primary TARGET \
    --primary-direction low \
    --output output/TARGET_tcga.json
```

**Arguments:**
- `--genes` — space-separated SYMBOL:ENTREZ pairs (e.g. `NMT2:10891 VIM:7431`)
- `--primary` — gene to sort cohorts by
- `--primary-direction` — `low` or `high`
- `--threshold` — z-score threshold (default 0.5)
- `--min-n` — minimum samples per cohort (default 30)

**Entrez IDs:** Look up at https://www.ncbi.nlm.nih.gov/gene or via cBioPortal gene search.

**Interpreting output:**
- `pct_low` — % samples with z < -threshold (marker-low / absent)
- `pct_high` — % samples with z > +threshold (marker-high / present)
- Cross-reference with co-markers (VIM, ZEB1) to distinguish mechanism-matched cohorts
  from false positives (e.g. AML: high marker-low% but hematopoietic lineage biology)

---

## Step 9 — BRCA PAM50 Subtype Stratification

Cross-reference cBioPortal PAM50 subtypes with PanCan RNA-seq z-scores.
Uses `brca_tcga_pub` for PAM50 annotation (522 samples) matched to PanCan RNA-seq.

```python
import urllib.request, json, statistics
from collections import defaultdict

API = 'https://www.cbioportal.org/api'
TARGET_ENTREZ = ENTREZ  # integer

# Get PAM50 subtypes
url = f'{API}/studies/brca_tcga_pub/clinical-data?attributeId=PAM50_SUBTYPE&projection=SUMMARY&pageSize=5000'
with urllib.request.urlopen(url) as r:
    subtype_map = {d['sampleId']: d['value'] for d in json.load(r)}

# Get RNA-seq z-scores from PanCan Atlas
profile = 'brca_tcga_pan_can_atlas_2018_rna_seq_v2_mrna_median_all_sample_Zscores'
payload = json.dumps({'entrezGeneIds': [TARGET_ENTREZ],
                      'sampleListId': 'brca_tcga_pan_can_atlas_2018_all'}).encode()
req = urllib.request.Request(f'{API}/molecular-profiles/{profile}/molecular-data/fetch',
                              data=payload, headers={'Content-Type': 'application/json'})
with urllib.request.urlopen(req, timeout=30) as r:
    mol_data = json.load(r)

vals_by_subtype = defaultdict(list)
gene_vals = {d['sampleId']: d['value'] for d in mol_data if d.get('value') is not None}
for sid, st in subtype_map.items():
    if sid in gene_vals:
        vals_by_subtype[st].append(gene_vals[sid])

for st in ['Basal-like','HER2-enriched','Luminal A','Luminal B','Normal-like']:
    vals = vals_by_subtype.get(st, [])
    if not vals: continue
    n = len(vals)
    print(f'{st}: n={n}  low={sum(v<-0.5 for v in vals)/n*100:.1f}%'
          f'  high={sum(v>0.5 for v in vals)/n*100:.1f}%'
          f'  med={statistics.median(vals):.3f}')
```

**Note:** `brca_tcga_pub` microarray z-scores are unusable for some genes (detection
floor). Always use PanCan RNA-seq (`_rna_seq_v2_mrna_median_all_sample_Zscores`) for
quantitative z-score analysis. PAM50 annotation from `brca_tcga_pub` matches ~99% of
PanCan samples by TCGA barcode.

---

## Interpreting the full picture

| Analysis | What it tells you |
|----------|------------------|
| XY Tier 2 (pearson_r) | Which expression features predict dependency; direction = biomarker type |
| XY p8_signed_wasserstein_shift | Positive = high marker → resistant (absence biomarker); negative = sensitive |
| XY p7_delta_aic | Threshold structure — sharp biomarker cut vs. continuous |
| XY p4_iqr_ratio | Variance modulation — heterogeneous response in marker-high cells |
| XX shape communities | Which co-expression modules underlie the XY biomarker communities |
| YY' positive correlates | Co-essential genes; same pathway/complex |
| YY' negative correlates | Sensitizer candidates — inhibit to create target dependency |
| YY' d_aic on anti-correlates | Threshold: sensitizer shifts cells past a regime boundary |
| TCGA prevalence | How many patients have the biomarker state in each cancer type |
| PAM50 stratification | Which breast cancer subtype is enriched for the biomarker |

### Patient selection logic

```
1. Identify primary biomarker from XY analysis
   - Absence biomarker (r < 0, wass > 0): target essential when marker is LOW
   - Presence biomarker (r > 0, wass < 0): target essential when marker is HIGH

2. Check TCGA prevalence: which cohorts have the highest biomarker-state %?
   Cross-reference with co-markers to confirm mechanistic plausibility.

3. Check PAM50 stratification for BRCA: which subtype dominates the signal?

4. Identify sensitizer combinations from YY' anti-correlates:
   - Druggable anti-correlates with high d_aic = combination leads
   - These expand the addressable patient population beyond the primary biomarker
```

---

## Output directory convention

```
output/TARGET_full/           # XY target mode outputs
output/TARGET_full/enrichr/   # Enrichr enrichment results
output/TARGET_full/joint_graph/  # Joint XX+XY graph
output/TARGET_full/xx_c1/     # XX Tier 2 on community C1 genes
output/TARGET_full/xx_c34/    # XX Tier 2 on community C3+C4 genes
output/TARGET_full/xx_combined_shape/  # Combined shape communities
output/TARGET_yy/             # YY' on expression-intersected chronos (n=277)
output/TARGET_yy_full/        # YY' on full chronos (n=1,538) -- preferred
output/TARGET_tcga.json       # TCGA pan-cancer prevalence
```

---

## NMT1 reference outputs

All NMT1 analysis outputs are in `output/NMT1_full/` and `output/NMT1_yy_full/`.
Full findings logged in `session_log_strategy_target_discovery.md`.
