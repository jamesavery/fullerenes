#!/usr/bin/env python3
"""
Analyze PipelineDiag bitfields from bench_epopt JSON output.

Cross-correlates diagnostic flags/counters with timing, convergence,
and quality metrics to identify which mechanisms predict slow convergence
and failures.

Usage:
  python3 diag_report.py bench_C60.json bench_C80.json ...
  python3 diag_report.py /tmp/diag_C*.json
"""
import json, sys, os
from collections import defaultdict
import math

# ── Flag definitions ──────────────────────────────────────────────────
FLAG_BITS = {
    0x01:      'reflect_cycling_step',
    0x02:      'hull_used_step',
    0x04:      'hull_cycling_step',
    0x08:      'cvx_fail_step',
    0x10:      'patch_cycling',
    0x20:      'stag_step',
    0x100:     'reflect_cycling_final',
    0x200:     'hull_used_final',
    0x400:     'hull_cycling_final',
    0x800:     'stag_final',
    0x1000:    'stag_constrained',
    0x2000:    'budget_constrained',
    0x10000:   'neg_curvature',
    0x20000:   'tr_boundary',
    0x40000:   'step_rejected',
    0x80000:   'cvx_rejected',
    0x100000:  'lbfgs_reset',
    0x1000000: 'has_f_ring',
}

ALL_FLAGS = [
    'reflect_cycling_step', 'hull_used_step', 'hull_cycling_step',
    'cvx_fail_step', 'patch_cycling', 'stag_step',
    'reflect_cycling_final', 'hull_used_final', 'hull_cycling_final',
    'stag_final', 'stag_constrained', 'budget_constrained',
    'neg_curvature', 'tr_boundary', 'step_rejected', 'cvx_rejected',
    'lbfgs_reset', 'has_f_ring',
]

COUNTER_NAMES = [
    'max_patch_rounds', 'max_reflect_rounds_step', 'final_reflect_rounds',
    'n_cvx_fail_steps', 'n_stag_steps',
]

ABBREVS = {
    'reflect_cycling_step': 'rfl_cyc_s', 'hull_used_step': 'hull_s',
    'hull_cycling_step': 'hcyc_s', 'cvx_fail_step': 'cvxf_s',
    'patch_cycling': 'pat_cyc', 'stag_step': 'stag_s',
    'reflect_cycling_final': 'rfl_cyc_f', 'hull_used_final': 'hull_f',
    'hull_cycling_final': 'hcyc_f', 'stag_final': 'stag_f',
    'stag_constrained': 'stag_c', 'budget_constrained': 'bud_c',
    'neg_curvature': 'neg_k', 'tr_boundary': 'tr_bnd',
    'step_rejected': 'rej', 'cvx_rejected': 'cvx_rej',
    'lbfgs_reset': 'lbfgs_r', 'has_f_ring': 'f_ring',
}

# ── Helpers ───────────────────────────────────────────────────────────
def decode_flags(hex_str):
    f = int(hex_str, 16)
    return {name for bit, name in FLAG_BITS.items() if f & bit}

def decode_counters(hex_str):
    c = int(hex_str, 16)
    return {
        'max_patch_rounds':        (c >> 0)  & 0xf,
        'max_reflect_rounds_step': (c >> 4)  & 0xf,
        'final_reflect_rounds':    (c >> 8)  & 0xf,
        'n_cvx_fail_steps':        (c >> 12) & 0xf,
        'n_stag_steps':            (c >> 16) & 0xf,
        'seed_type':               (c >> 24) & 0x3,
    }

def median(vals):
    s = sorted(vals)
    n = len(s)
    if n == 0: return 0
    return s[n // 2]

def mean(vals):
    return sum(vals) / len(vals) if vals else 0

def pct(vals, p):
    s = sorted(vals)
    n = len(s)
    if n == 0: return 0
    return s[min(int(n * p / 100), n - 1)]

def fmt_pct(count, total):
    """Format count/total as percentage or raw count."""
    if count == 0: return '·'.rjust(9)
    if count == total: return 'ALL'.rjust(9)
    p = 100.0 * count / total
    if p < 1: return f'{count}'.rjust(9)
    return f'{p:.1f}%'.rjust(9)

# ── Load data ─────────────────────────────────────────────────────────
files = sys.argv[1:]
if not files:
    print("Usage: diag_report.py bench_C60.json bench_C80.json ...")
    sys.exit(1)

by_size = defaultdict(list)  # N -> list of (isomer_dict, flags_set, counters_dict)

for path in files:
    try:
        d = json.load(open(path))
    except Exception as e:
        print(f"Warning: cannot read {path}: {e}", file=sys.stderr)
        continue
    N = d['N']
    for p in d['per_isomer']:
        if 'diag_flags' not in p:
            continue
        flags = decode_flags(p['diag_flags'])
        counters = decode_counters(p['diag_counters'])
        by_size[N].append((p, flags, counters))

total = sum(len(v) for v in by_size.values())
if total == 0:
    print("No isomers with diag data found.")
    sys.exit(1)
print(f"Total isomers: {total}")
print()

# ══════════════════════════════════════════════════════════════════════
# 1. Per-size flag prevalence table
# ══════════════════════════════════════════════════════════════════════
print("=" * 100)
print("1. FLAG PREVALENCE BY SIZE")
print("=" * 100)
# Only show flags that appear at least once
active_flags = [fn for fn in ALL_FLAGS
                if any(fn in f for items in by_size.values() for _, f, _ in items)]

print(f"{'Size':>5} {'n':>6} {'!cv':>4} |", end='')
for fn in active_flags:
    print(f' {ABBREVS.get(fn, fn[:8]):>9}', end='')
print()
print("-" * (18 + 10 * len(active_flags)))

for N in sorted(by_size.keys()):
    items = by_size[N]
    n_total = len(items)
    n_nc = sum(1 for p, f, c in items if not p['converged'])
    print(f"C{N:>3} {n_total:>6} {n_nc:>4} |", end='')
    for fn in active_flags:
        count = sum(1 for p, f, c in items if fn in f)
        print(f' {fmt_pct(count, n_total)}', end='')
    print()

# ══════════════════════════════════════════════════════════════════════
# 2. Flag-timing correlation: median ms with vs without each flag
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 100)
print("2. FLAG-TIMING CORRELATION (median ms with/without flag, sorted by ratio)")
print("=" * 100)
print(f"{'Size':>5} | {'Flag':>25} | {'n_with':>6} {'med_with':>10} {'p95_with':>10}"
      f" | {'n_without':>9} {'med_without':>11} {'ratio':>7}")
print("-" * 100)

for N in sorted(by_size.keys()):
    items = by_size[N]
    rows = []
    for fn in active_flags:
        wf  = [p['ms'] for p, f, c in items if fn in f]
        wof = [p['ms'] for p, f, c in items if fn not in f]
        if len(wf) < 1 or len(wof) < 1:
            continue
        ratio = median(wf) / max(1e-9, median(wof))
        # Skip universal flags with no timing difference
        if ratio < 1.03 and len(wf) == len(items):
            continue
        rows.append((fn, wf, wof, ratio))
    rows.sort(key=lambda x: -x[3])
    for fn, wf, wof, ratio in rows:
        print(f"C{N:>3} | {fn:>25} | {len(wf):>6} {median(wf):>10.0f} {pct(wf, 95):>10.0f}"
              f" | {len(wof):>9} {median(wof):>11.0f} {ratio:>6.2f}x")

# ══════════════════════════════════════════════════════════════════════
# 3. Counter-timing correlation
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 100)
print("3. COUNTER-TIMING CORRELATION (median ms by counter value)")
print("=" * 100)

for N in sorted(by_size.keys()):
    items = by_size[N]
    for cn in COUNTER_NAMES:
        by_val = defaultdict(list)
        for p, f, c in items:
            by_val[c[cn]].append(p['ms'])
        vals_present = sorted(by_val.keys())
        if len(vals_present) <= 1 or max(vals_present) == 0:
            continue
        print(f"C{N:>3} {cn:>30}:", end='')
        for v in vals_present:
            ms_list = by_val[v]
            print(f"  {v}->{median(ms_list):.0f}ms(n={len(ms_list)})", end='')
        print()

# ══════════════════════════════════════════════════════════════════════
# 4. Slow-isomer enrichment (top 5% by ms)
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 100)
print("4. SLOW-ISOMER ENRICHMENT (flags enriched in top 5% slowest)")
print("=" * 100)

for N in sorted(by_size.keys()):
    items = by_size[N]
    if len(items) < 20:
        continue
    ms_vals = sorted([p['ms'] for p, f, c in items])
    threshold = ms_vals[int(len(ms_vals) * 0.95)]
    slow = [(p, f, c) for p, f, c in items if p['ms'] >= threshold]
    fast = [(p, f, c) for p, f, c in items if p['ms'] < threshold]
    if not slow:
        continue

    print(f"\nC{N}: top 5% = {len(slow)} isomers, ms >= {threshold:.0f}")
    rows = []
    for fn in active_flags:
        n_slow = sum(1 for p, f, c in slow if fn in f)
        n_fast = sum(1 for p, f, c in fast if fn in f)
        pct_slow = 100.0 * n_slow / len(slow)
        pct_fast = 100.0 * n_fast / len(fast) if fast else 0
        enrichment = pct_slow / max(0.1, pct_fast) if pct_fast > 0 else (999 if n_slow > 0 else 0)
        if enrichment > 1.5 or (n_slow > 0 and n_fast == 0):
            rows.append((fn, pct_slow, pct_fast, enrichment))
    rows.sort(key=lambda x: -x[3])
    for fn, ps, pf, enr in rows:
        print(f"  {fn:>30}: {ps:5.1f}% slow vs {pf:5.1f}% fast  ({enr:.1f}x)")

# ══════════════════════════════════════════════════════════════════════
# 5. Worst outliers per size
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 100)
print("5. WORST OUTLIERS BY ms (top 5 per size)")
print("=" * 100)

for N in sorted(by_size.keys()):
    items = by_size[N]
    items_sorted = sorted(items, key=lambda x: -x[0]['ms'])
    print(f"\nC{N}:")
    for p, flags, counters in items_sorted[:5]:
        flags_str = ', '.join(sorted(flags - {'patch_cycling'})) or '(clean)'
        ct_parts = []
        for k in COUNTER_NAMES:
            v = counters[k]
            if v > 0:
                ct_parts.append(f'{k}={v}')
        ct_str = ', '.join(ct_parts) or '-'
        print(f"  idx={p['idx']:>7} {p['ms']:>8.0f}ms {p.get('exit','?'):>17}"
              f"  ang={p['ang_relerr_max']:.1e} gmax={p['gmax_L']:.1e}"
              f"  [{flags_str}]  ({ct_str})")

# ══════════════════════════════════════════════════════════════════════
# 6. Seed-type breakdown
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 100)
print("6. SEED-TYPE BREAKDOWN")
print("=" * 100)

seed_names = {0: 'C20', 1: 'C28', 2: 'C30'}
for N in sorted(by_size.keys()):
    items = by_size[N]
    by_seed = defaultdict(list)
    for p, f, c in items:
        by_seed[c['seed_type']].append(p)
    if len(by_seed) <= 1:
        continue
    parts = []
    for st in sorted(by_seed.keys()):
        ps = by_seed[st]
        n = len(ps)
        med = median([p['ms'] for p in ps])
        nc = sum(1 for p in ps if not p['converged'])
        parts.append(f"{seed_names.get(st, f'?{st}')}:{n}({med:.0f}ms, !cv={nc})")
    print(f"C{N:>3}: {', '.join(parts)}")

# ══════════════════════════════════════════════════════════════════════
# 7. Work efficiency: work per vertex^2 by flag
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 100)
print("7. WORK EFFICIENCY (median work/Nv^2 with vs without flag)")
print("=" * 100)

for N in sorted(by_size.keys()):
    items = by_size[N]
    Nv = N // 2 + 2
    rows = []
    for fn in active_flags:
        wf  = [p['work'] / (Nv * Nv) for p, f, c in items if fn in f and p.get('work', 0) > 0]
        wof = [p['work'] / (Nv * Nv) for p, f, c in items if fn not in f and p.get('work', 0) > 0]
        if len(wf) < 1 or len(wof) < 1:
            continue
        ratio = median(wf) / max(1e-9, median(wof))
        if abs(ratio - 1.0) < 0.03 and len(wf) == len(items):
            continue
        rows.append((fn, len(wf), median(wf), len(wof), median(wof), ratio))
    rows.sort(key=lambda x: -x[5])
    for fn, nw, mw, nwo, mwo, ratio in rows[:8]:
        print(f"C{N:>3} {fn:>25}: {mw:>7.1f} (n={nw:>5}) vs {mwo:>7.1f} (n={nwo:>5})  {ratio:.2f}x")

# ══════════════════════════════════════════════════════════════════════
# 8. Quality degradation: angle error by flag
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 100)
print("8. QUALITY: median ang_relerr_max with vs without flag")
print("=" * 100)

for N in sorted(by_size.keys()):
    items = by_size[N]
    rows = []
    for fn in active_flags:
        wf  = [p['ang_relerr_max'] for p, f, c in items if fn in f]
        wof = [p['ang_relerr_max'] for p, f, c in items if fn not in f]
        if len(wf) < 1 or len(wof) < 1:
            continue
        mw, mwo = median(wf), median(wof)
        if mwo < 1e-20:
            continue
        ratio = mw / mwo
        if abs(ratio - 1.0) < 0.1 and len(wf) == len(items):
            continue
        rows.append((fn, len(wf), mw, len(wof), mwo, ratio))
    rows.sort(key=lambda x: -x[5])
    for fn, nw, mw, nwo, mwo, ratio in rows[:6]:
        print(f"C{N:>3} {fn:>25}: {mw:.1e} (n={nw:>5}) vs {mwo:.1e} (n={nwo:>5})  {ratio:.1f}x")

# ══════════════════════════════════════════════════════════════════════
# 9. Flag combination analysis (most common flag sets among outliers)
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 100)
print("9. FLAG COMBINATIONS (most common among top-10% slowest)")
print("=" * 100)

for N in sorted(by_size.keys()):
    items = by_size[N]
    if len(items) < 50:
        continue
    ms_vals = sorted([p['ms'] for p, f, c in items])
    threshold = ms_vals[int(len(ms_vals) * 0.90)]
    slow = [(p, f, c) for p, f, c in items if p['ms'] >= threshold]

    combo_counts = defaultdict(int)
    for p, flags, counters in slow:
        # Remove patch_cycling (universal) for cleaner combos
        key = frozenset(flags - {'patch_cycling'})
        combo_counts[key] += 1

    ranked = sorted(combo_counts.items(), key=lambda x: -x[1])
    print(f"\nC{N}: top 10% = {len(slow)} isomers")
    for combo, count in ranked[:5]:
        pct_val = 100.0 * count / len(slow)
        combo_str = ', '.join(sorted(combo)) if combo else '(clean)'
        print(f"  {count:>4} ({pct_val:4.1f}%): {combo_str}")

# ══════════════════════════════════════════════════════════════════════
# 10. Scaling trends: how flag prevalence grows with N
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 100)
print("10. SCALING: flag prevalence vs N")
print("=" * 100)

# Show flags whose prevalence changes significantly with N
sizes = sorted(by_size.keys())
if len(sizes) >= 3:
    for fn in active_flags:
        rates = []
        for N in sizes:
            items = by_size[N]
            count = sum(1 for p, f, c in items if fn in f)
            rates.append(100.0 * count / len(items))
        # Skip if constant (all 0 or all 100) or nearly so
        if max(rates) - min(rates) < 5:
            continue
        parts = [f"C{N}:{r:.0f}%" for N, r in zip(sizes, rates)]
        print(f"  {fn:>30}: {', '.join(parts)}")
