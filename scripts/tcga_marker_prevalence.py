"""
TCGA pan-cancer marker prevalence via cBioPortal API.

For a set of genes, fetches pan-cancer z-score expression across all 30 TCGA
PanCancer Atlas 2018 cohorts and computes:
  - % samples with z < threshold (default -0.5 = "low")
  - % samples with z > threshold (default +0.5 = "high")
  - median z-score per cohort

Useful for patient selection strategy development: which cancer types have the
highest prevalence of a given marker state?

Usage:
    # NMT2-low + VIM/ZEB1 cross-reference
    python scripts/tcga_marker_prevalence.py \
        --genes NMT2:10891 VIM:7431 ZEB1:6935 \
        --primary NMT2 --primary-direction low \
        --output output/nmt2_vim_zeb1_tcga.json

    # Any target: single gene, high direction
    python scripts/tcga_marker_prevalence.py \
        --genes MTOR:2475 \
        --primary MTOR --primary-direction high \
        --output output/mtor_tcga.json
"""

import argparse
import json
import statistics
import time
import urllib.request
from pathlib import Path


TCGA_PANCAN_SUFFIX = "_tcga_pan_can_atlas_2018"
ZSCORE_PROFILE_SUFFIX = "_rna_seq_v2_mrna_median_all_sample_Zscores"
API_BASE = "https://www.cbioportal.org/api"


def fetch_zscores(molecular_profile_id: str, sample_list_id: str, entrez_id: int) -> list[float]:
    url = f"{API_BASE}/molecular-profiles/{molecular_profile_id}/molecular-data/fetch"
    payload = json.dumps({"entrezGeneIds": [entrez_id], "sampleListId": sample_list_id}).encode()
    req = urllib.request.Request(url, data=payload, headers={"Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=20) as r:
            return [d["value"] for d in json.load(r) if d.get("value") is not None]
    except Exception as e:
        return []


def get_tcga_studies() -> list[dict]:
    url = f"{API_BASE}/studies?pageSize=500&projection=SUMMARY"
    with urllib.request.urlopen(url) as r:
        studies = json.load(r)
    return [s for s in studies if "tcga_pan_can_atlas" in s["studyId"]]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genes", nargs="+", required=True,
                        help="SYMBOL:ENTREZ pairs, e.g. NMT2:10891 VIM:7431")
    parser.add_argument("--primary", required=True,
                        help="Primary gene symbol to sort by")
    parser.add_argument("--primary-direction", choices=["low", "high"], default="low",
                        help="Sort by %% low (z < -0.5) or %% high (z > 0.5)")
    parser.add_argument("--threshold", type=float, default=0.5,
                        help="Absolute z-score threshold (default 0.5)")
    parser.add_argument("--min-n", type=int, default=30,
                        help="Minimum samples per cohort (default 30)")
    parser.add_argument("--output", required=True,
                        help="Output JSON path")
    parser.add_argument("--delay", type=float, default=0.15,
                        help="Seconds between API calls (default 0.15)")
    args = parser.parse_args()

    genes = {}
    for g in args.genes:
        sym, eid = g.split(":")
        genes[sym] = int(eid)

    print(f"Genes: {genes}")
    print(f"Fetching TCGA PanCan studies...")
    studies = get_tcga_studies()
    print(f"  {len(studies)} studies found")

    results = []
    for s in studies:
        sid = s["studyId"]
        name = s["name"].split("(")[0].strip()
        zpid = sid + ZSCORE_PROFILE_SUFFIX
        slist = sid + "_all"

        row = {"cancer": name, "study": sid}
        for sym, entrez in genes.items():
            vals = fetch_zscores(zpid, slist, entrez)
            time.sleep(args.delay)
            if len(vals) < args.min_n:
                row[f"{sym}_n"] = len(vals)
                row[f"{sym}_pct_low"] = None
                row[f"{sym}_pct_high"] = None
                row[f"{sym}_median_z"] = None
            else:
                n = len(vals)
                row[f"{sym}_n"] = n
                row[f"{sym}_pct_low"]  = round(sum(1 for v in vals if v < -args.threshold) / n * 100, 1)
                row[f"{sym}_pct_high"] = round(sum(1 for v in vals if v >  args.threshold) / n * 100, 1)
                row[f"{sym}_median_z"] = round(statistics.median(vals), 3)

        if row.get(f"{args.primary}_n", 0) >= args.min_n:
            results.append(row)
            parts = [f"{sym}={'%.1f' % row[f'{sym}_pct_low']}%low {'%.3f' % row[f'{sym}_median_z']}z"
                     for sym in genes if row.get(f'{sym}_n', 0) >= args.min_n]
            print(f"  {sid}: {' | '.join(parts)}")

    sort_key = f"{args.primary}_pct_{args.primary_direction}"
    results.sort(key=lambda x: -(x.get(sort_key) or 0))

    # Print ranked table
    col_w = 44
    headers = ["Cancer".ljust(col_w)] + [f"{sym} {'low%' if args.primary_direction=='low' else 'high%':>8} med_z" for sym in genes]
    print("\n=== RANKED BY " + args.primary + "-" + args.primary_direction.upper() + " ===")
    print("  ".join(headers))
    print("-" * (col_w + 20 * len(genes)))
    for r in results:
        row_str = r["cancer"].ljust(col_w)
        for sym in genes:
            pct = r.get(f"{sym}_pct_{args.primary_direction}")
            mz  = r.get(f"{sym}_median_z")
            pct_s = f"{pct:>6.1f}%" if pct is not None else "    n/a"
            mz_s  = f"{mz:>7.3f}"  if mz  is not None else "    n/a"
            row_str += f"  {pct_s} {mz_s}"
        print(row_str)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved -> {out}  ({len(results)} cohorts)")


if __name__ == "__main__":
    main()
