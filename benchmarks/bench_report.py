#!/usr/bin/env python3
"""Format bench_epopt JSON output as a convergence-report-style table."""
import json, sys, math

def fmt_sci(x, digits=1):
    if x == 0: return "0"
    if math.isnan(x) or math.isinf(x): return str(x)
    exp = math.floor(math.log10(abs(x)))
    mantissa = x / 10**exp
    if digits == 1:
        return f"{mantissa:.0f}e{exp:+d}" if exp != 0 else f"{mantissa:.0f}"
    return f"{mantissa:.{digits-1}f}e{exp:+d}"

def fmt_pct(x):
    """Format as percentage with appropriate precision."""
    if x < 0.01: return f"{x*100:.2f}%"
    if x < 0.1: return f"{x*100:.1f}%"
    return f"{x*100:.0f}%"

def get_ang_relerr(iso):
    """Get angle relative error, computing from ang_min/ang_max if not stored."""
    if "ang_relerr_max" in iso:
        return iso["ang_relerr_max"]
    return max(60.0 - iso["ang_min"], iso["ang_max"] - 60.0) / 60.0

def get_edge_relerr(iso):
    """Get edge relative error if stored, else None."""
    return iso.get("edge_relerr_max", None)

def report(path):
    with open(path) as f:
        d = json.load(f)

    N = d["N"]
    M = d["M"]
    dual_N = d["dual_N"]
    total_enum = d.get("total_enumerated", M)
    s = d["summary"]
    total_s = d["total_ms"] / 1000
    isomers = d["per_isomer"]

    # Compute per-isomer derived metrics
    for iso in isomers:
        iso["_ang_relerr"] = get_ang_relerr(iso)
        iso["_edge_relerr"] = get_edge_relerr(iso)

    # Compute summary metrics not in JSON
    ang_relerr_worst = max(iso["_ang_relerr"] for iso in isomers)
    edge_relerr_worst = None
    if isomers[0]["_edge_relerr"] is not None:
        edge_relerr_worst = max(iso["_edge_relerr"] for iso in isomers)

    n_not_conv = s.get("n_not_converged", sum(1 for iso in isomers if not iso["converged"]))
    n_bad_edge = s.get("n_bad_edge", None)
    n_bad_ang = s.get("n_bad_ang", sum(1 for iso in isomers if iso["_ang_relerr"] > 0.01))

    sample_note = f"{M}/{total_enum}" if total_enum > M else f"all {M}"
    print(f"## C{N} (dual N={dual_N}, {sample_note} isomers, {total_s:.1f}s)")
    print()
    print(f"| Metric | Mean | Worst | Bad (>1%) |")
    print(f"|---|---|---|---|")
    print(f"| ms/isomer | {s['ms_mean']:.1f} +/- {s['ms_std']:.1f} | - | - |")
    print(f"| iters/vertex | {s['iters_per_vertex']:.1f} | - | - |")
    print(f"| edge CV | {fmt_sci(s['edge_cv_mean'])} | {fmt_sci(s['edge_cv_max'])} | - |")
    if edge_relerr_worst is not None:
        print(f"| edge max relerr | - | {fmt_pct(edge_relerr_worst)} | {n_bad_edge}/{M} |")
    print(f"| ang CV | {fmt_sci(s['ang_std_mean'] / 60)} | {fmt_sci(s['ang_std_max'] / 60)} | - |")
    print(f"| ang max relerr | - | {fmt_pct(ang_relerr_worst)} | {n_bad_ang}/{M} |")
    print(f"| ang range | - | [{s['ang_min_worst']:.1f}, {s['ang_max_worst']:.1f}] | - |")
    print(f"| h_min | - | {s['h_min_worst']:+.4f} | - |")
    print(f"| concave | - | - | {s['concave_isomers']}/{M} |")
    print(f"| gmax*L | {fmt_sci(s['gmax_L_mean'])} | {fmt_sci(s['gmax_L_max'])} | {n_not_conv}/{M} not conv |")

    # Worst 5 by gmax*L
    by_gmax = sorted(isomers, key=lambda x: -x["gmax_L"])[:5]
    print()
    print(f"Worst 5 by gmax*L:")
    hdr =  "| idx | ms | iters | edge_cv | edge_re | ang_re | h_min | conc | gmax*L | seed | steps |"
    sep =  "|---|---|---|---|---|---|---|---|---|---|---|"
    print(hdr)
    print(sep)
    for r in by_gmax:
        er = fmt_pct(r["_edge_relerr"]) if r["_edge_relerr"] is not None else "-"
        print(f"| {r['idx']} | {r['ms']:.0f} | {r['iters']} | {fmt_sci(r['edge_cv'])} "
              f"| {er} | {fmt_pct(r['_ang_relerr'])} "
              f"| {r['h_min']:+.3f} | {r['n_concave']} | {fmt_sci(r['gmax_L'])} "
              f"| {r['seed']} | {r['n_steps']} |")

    # Worst 5 by angle relerr
    by_ang = sorted(isomers, key=lambda x: -x["_ang_relerr"])[:5]
    print()
    print(f"Worst 5 by angle relerr:")
    print(hdr)
    print(sep)
    for r in by_ang:
        er = fmt_pct(r["_edge_relerr"]) if r["_edge_relerr"] is not None else "-"
        print(f"| {r['idx']} | {r['ms']:.0f} | {r['iters']} | {fmt_sci(r['edge_cv'])} "
              f"| {er} | {fmt_pct(r['_ang_relerr'])} "
              f"| {r['h_min']:+.3f} | {r['n_concave']} | {fmt_sci(r['gmax_L'])} "
              f"| {r['seed']} | {r['n_steps']} |")
    print()

def summary_table(paths):
    """Print a compact comparison table across multiple sizes."""
    rows = []
    for path in paths:
        with open(path) as f:
            d = json.load(f)
        s = d["summary"]
        total_enum = d.get("total_enumerated", d["M"])
        isomers = d["per_isomer"]

        # Compute derived metrics
        ang_re_max = max(max(60.0 - iso["ang_min"], iso["ang_max"] - 60.0) / 60.0 for iso in isomers)
        edge_re_max = s.get("edge_relerr_max", None)
        n_not_conv = s.get("n_not_converged", sum(1 for iso in isomers if not iso["converged"]))
        n_bad_ang = s.get("n_bad_ang", sum(1 for iso in isomers if
                          max(60.0 - iso["ang_min"], iso["ang_max"] - 60.0) / 60.0 > 0.01))
        n_bad_edge = s.get("n_bad_edge", None)

        rows.append({
            "N": d["N"], "dual_N": d["dual_N"], "M": d["M"],
            "total": total_enum,
            "ang_re_max": ang_re_max, "edge_re_max": edge_re_max,
            "n_not_conv": n_not_conv, "n_bad_ang": n_bad_ang, "n_bad_edge": n_bad_edge,
            **s
        })

    print("## Summary across sizes")
    print()
    print("| N | dual | M/total | ms/iso | it/vtx | edge_cv | edge_re | ang_cv | ang_re | h_min | conc | !conv | !edge | !ang |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        sample = f"{r['M']}/{r['total']}" if r['total'] > r['M'] else f"{r['M']}"
        edge_re = fmt_pct(r['edge_re_max']) if r['edge_re_max'] is not None else "-"
        n_be = str(r['n_bad_edge']) if r['n_bad_edge'] is not None else "-"
        ang_cv = fmt_sci(r['ang_std_mean'] / 60)
        print(f"| C{r['N']} | {r['dual_N']} | {sample} "
              f"| {r['ms_mean']:.0f} | {r['iters_per_vertex']:.0f} "
              f"| {fmt_sci(r['edge_cv_max'])} "
              f"| {edge_re} "
              f"| {ang_cv} "
              f"| {fmt_pct(r['ang_re_max'])} "
              f"| {r['h_min_worst']:+.3f} "
              f"| {r['concave_isomers']} "
              f"| {r['n_not_conv']} "
              f"| {n_be} "
              f"| {r['n_bad_ang']} |")
    print()

if __name__ == "__main__":
    paths = sys.argv[1:]
    if len(paths) > 1:
        summary_table(paths)
    for path in paths:
        report(path)
