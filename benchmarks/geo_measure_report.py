#!/usr/bin/env python3
"""Tables for the .geo format decisions from geo_measure JSON lines.

Usage: geo_measure_report.py C60.jsonl C80.jsonl ...

Per fullerene size N and embedding (dual: 12 cones; cubic: 20..60 cones), prints the
graph statistics and, for fixed-point widths w, bytes per record against the
worst-case position error (Euclidean, in edge lengths) for these coordinate framings:

  bound  one scale per file from a proven extent bound (safe for appending):
         |x_u - x_v| <= d_G(u,v) since every edge is an intrinsic unit segment, so a
         bounding-box-centred record has half-extent <= diam/2, with diam <= N/5 + 1
         for the cubic graph (Andova et al. 2012) and <= N/5 + 2 for the dual
         (d_T <= d_G + 1: consecutive faces along a cubic path share a vertex)
  batch  one scale per file from the largest half-extent in this sample
  rec    one f32 scale per record (+4 B); the caller centres the bounding box
  axes   three f32 scales per record (+12 B); the caller centres the bounding box,
         in the given or the principal frame, whichever is better
Storing the offset instead of centring would cost another 12 B per record.

Error bound: sqrt(3)*s/2 for one scale s, sqrt(sum s_k^2)/2 for per-axis scales,
where s = half-extent / (2^(w-1) - 1).  Records are rounded up to whole bytes and
then to a multiple of 4.  Graph: DCEL twin matching with fixed capacities --
n_cap * deg_bits + (A_cap/2) * bit_length(A_cap - 1), plus the per-record cone
count for the cubic embedding; "var" is the variable-length alternative with a
64-bit index entry per record.
"""
import json
import math
import statistics
import sys
from collections import defaultdict

WIDTHS = (8, 10, 12, 14, 16)


def bits_for(values):            # bits to store a value in [0, values)
    return (values - 1).bit_length()


def record_bytes(bits):
    return 4 * math.ceil(math.ceil(bits / 8) / 4)


def graph_layout(kind, recs):
    dmin = min(r["dmin"] for r in recs)
    dmax = max(r["dmax"] for r in recs)
    deg_bits = bits_for(dmax - dmin + 1)
    n_cap = 12 if kind == "dual" else 60
    A_cap = min(n_cap * dmax, 6 * n_cap - 12)
    cap_bits = n_cap * deg_bits + (A_cap // 2) * bits_for(A_cap)
    if kind == "cubic":
        cap_bits += bits_for(n_cap + 1)
    return dmin, dmax, deg_bits, n_cap, cap_bits


def quantiles(xs):
    xs = sorted(xs)
    return xs[0], statistics.median(xs), xs[-1]


def report(N, kind, recs, errs):
    print(f"\n### C{N} {kind}: {len(recs)} ok, {len(errs)} failed")
    for msg, count in sorted(defaultdict_count(errs).items(), key=lambda kv: -kv[1])[:5]:
        print(f"  failure x{count}: {msg[:150]}")
    if not recs:
        return

    dmin, dmax, deg_bits, n_cap, cap_bits = graph_layout(kind, recs)
    n_lo, n_med, n_hi = quantiles([r["n"] for r in recs])
    nonsimple = sum(1 for r in recs if r["loops"] or r["multi"])
    diam = max(r["diam"] for r in recs)
    bound_diam = N / 5 + (1 if kind == "cubic" else 2)
    print(f"  cones n: min {n_lo}, median {n_med}, max {n_hi};  degrees {dmin}..{dmax} "
          f"({deg_bits} bits);  non-simplicial iDT: {nonsimple}")
    print(f"  graph diameter: max {diam} (bound {bound_diam:g});  "
          f"graph block at capacity n_cap={n_cap}: {cap_bits} bits")
    if kind == "dual":
        lo, med, hi = quantiles([r["tbar_arcs"] for r in recs])
        print(f"  T-bar arcs: min {lo}, median {med}, max {hi} (capacity 60)")

    q = {w: 2 ** (w - 1) - 1 for w in WIDTHS}
    R_bound = bound_diam / 2
    R_batch = max(r["hb"][0] for r in recs)

    def axes_err(r, w):
        return min(math.sqrt(sum((h / q[w]) ** 2 for h in r[f])) / 2 for f in ("hb", "hp"))

    print(f"\n  | w | layout | bytes/rec | bound max err | batch max err "
          f"| rec max / median err | axes max / median err |")
    print(f"  |---|---|---|---|---|---|---|")
    for w in WIDTHS:
        coord_cap = 3 * n_cap * w
        layouts = [("cap", lambda r: coord_cap + cap_bits)]
        if kind == "cubic":
            layouts.append(("var", lambda r: 3 * r["n"] * w + bits_for(n_cap + 1) + r["n"] * deg_bits
                            + (r["A"] // 2) * bits_for(r["A"]) + 64))
        rec_err = [math.sqrt(3) * r["hb"][0] / q[w] / 2 for r in recs]
        axes = [axes_err(r, w) for r in recs]
        for name, bits in layouts:
            base = statistics.mean(record_bytes(bits(r)) for r in recs)
            rec = statistics.mean(record_bytes(bits(r) + 32) for r in recs)
            ax = statistics.mean(record_bytes(bits(r) + 96) for r in recs)
            print(f"  | {w} | {name} | {base:.0f} / {rec:.0f} / {ax:.0f} "
                  f"| {math.sqrt(3) * R_bound / q[w] / 2:.2e} "
                  f"| {math.sqrt(3) * R_batch / q[w] / 2:.2e} "
                  f"| {max(rec_err):.2e} / {statistics.median(rec_err):.2e} "
                  f"| {max(axes):.2e} / {statistics.median(axes):.2e} |")
    print("  (bytes/rec: file scale / + per-record scale / + per-axis scales)")


def defaultdict_count(msgs):
    d = defaultdict(int)
    for m in msgs:
        d[m] += 1
    return d


def main(paths):
    data = defaultdict(lambda: defaultdict(lambda: ([], [])))
    for path in paths:
        with open(path) as f:
            for line in f:
                r = json.loads(line)
                for kind in ("dual", "cubic"):
                    recs, errs = data[r["N"]][kind]
                    if kind in r:
                        recs.append(r[kind])
                    else:
                        errs.append(r[kind + "_err"])
    for N in sorted(data):
        for kind in ("dual", "cubic"):
            report(N, kind, *data[N][kind])


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1:])
