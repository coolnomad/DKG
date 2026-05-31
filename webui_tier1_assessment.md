Let me actually measure it on representative data.Single-core sandbox (the user's machine almost certainly has more). I'll benchmark here and we can scale.Interesting result: **`scipy.stats.t.sf` is 87% of the runtime, not the matmul**. Let me check whether p-values are even needed at full resolution (since the filter is on |r|, not p):Concrete answer based on the actual shapes:

**On 1 core (this sandbox), current code: ~2 minutes.**

But the bottleneck is not what you'd expect. Here's the breakdown:

| stage | time | share |
|---|---|---|
| matmul (Pearson + Spearman blocks) | 11.4 s | 9% |
| `scipy.stats.t.sf` on full blocks | **107 s** | **87%** |
| rank-transform | 1.8 s | 1% |
| centering + everything else | ~1 s | <1% |
| **total** | **~122 s** | |

87% of the time is computing p-values for pairs you're about to throw away. The filter is on `|r| >= threshold`, but `_pearson_block` evaluates `t.sf` for all 500 × 11,744 entries on every chunk — ~480 million survival-function calls in total, only ~150 of which actually survive at threshold 0.3.

**Easy fix: defer p-value computation to surviving pairs.** Compute `r` for the whole block, threshold, then run `t.sf` only on the survivors. Measured here:

| threshold | survivors | runtime |
|---|---|---|
| 0.30 | 156 | **13.1 s** |
| 0.20 | 285k | 13.3 s |
| 0.15 | 4.0M | 15.4 s |
| 0.10 | 29.0M | 32.0 s |

So a ~9× speedup on this machine, with t.sf cost now negligible unless your threshold lets through millions of pairs.

The patch is small — in `_pearson_block`, return `r` only; do the t→p step inside the `_screen_no_na` loop after `np.where(mask)`:

```python
# After: ci, cj = np.where(mask)
pr_s, sr_s = pr[ci, cj], sr[ci, cj]
tp = pr_s * np.sqrt((n - 2) / np.maximum(1.0 - pr_s**2, 1e-300))
ts = sr_s * np.sqrt((n - 2) / np.maximum(1.0 - sr_s**2, 1e-300))
pp = 2.0 * scipy.stats.t.sf(np.abs(tp), df=n - 2)
sp = 2.0 * scipy.stats.t.sf(np.abs(ts), df=n - 2)
```

**Scaling to your actual machine:** This sandbox is 1 core. The matmul is the one piece that scales nearly linearly via OpenBLAS — on 8 cores you'd see ~1.5 s instead of 11.4 s. Rank-transform doesn't parallelize as cleanly in scipy and will likely stay close to its serial cost. t.sf, after the fix, is irrelevant.

Rough projection for your run:

| | 1 core (measured) | 4 cores | 8+ cores |
|---|---|---|---|
| current code | ~120 s | ~50–60 s | ~30–40 s |
| with deferred t.sf | ~13 s | ~5–7 s | **~3–5 s** |

If you're on a recent laptop/workstation with the fix applied, expect well under 10 seconds. If you're seeing minutes in practice, the t.sf overhead is almost certainly the reason.

One nuance for the deferred-p-value version: at the very loose-threshold end (`thr=0.10` with 29M survivors), the t.sf cost comes back. That's still the right tradeoff — you only pay for p-values you actually want — but it's worth knowing the runtime is now threshold-dependent in a way the original code wasn't.