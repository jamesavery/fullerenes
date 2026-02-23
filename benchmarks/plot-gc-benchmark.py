#!/usr/bin/env python3
"""Plot GC transform benchmark results."""

import json
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) < 2:
    print(f"Usage: {sys.argv[0]} <benchmark.json>", file=sys.stderr)
    sys.exit(1)

with open(sys.argv[1]) as f:
    data = json.load(f)

out_dir = os.path.dirname(os.path.abspath(sys.argv[1]))

halma = data["halma"]
chiral = data["chiral"]

h_N = np.array([d["N_carbon"] for d in halma])
h_t = np.array([d["time_us"] for d in halma])
c_N = np.array([d["N_carbon"] for d in chiral])
c_t = np.array([d["time_us"] for d in chiral])

# Convert to milliseconds for readability
h_t_ms = h_t / 1000.0
c_t_ms = c_t / 1000.0

# --- Plot 1: Time per graph ---
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(h_N, h_t_ms, 'o-', label='Halma (l=0)', markersize=4)
ax.plot(c_N, c_t_ms, 's', label='Chiral (l>0)', markersize=3, alpha=0.6)
ax.set_xlabel('N (carbon atoms)')
ax.set_ylabel('Time (ms)')
ax.set_title('GC Transform Time — C28 base')
ax.legend()
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(out_dir, 'gc-benchmark-time.png'), dpi=150)
print(f"Saved {os.path.join(out_dir, 'gc-benchmark-time.png')}")

# --- Plot 2: Time per vertex ---
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(h_N, h_t / h_N, 'o-', label='Halma (l=0)', markersize=4)
ax.plot(c_N, c_t / c_N, 's', label='Chiral (l>0)', markersize=3, alpha=0.6)
ax.set_xlabel('N (carbon atoms)')
ax.set_ylabel('Time per vertex (μs)')
ax.set_title('GC Transform Time per Vertex — C28 base')
ax.legend()
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(out_dir, 'gc-benchmark-per-vertex.png'), dpi=150)
print(f"Saved {os.path.join(out_dir, 'gc-benchmark-per-vertex.png')}")
