#!/usr/bin/env python3
"""Plot spiral benchmark results: computation and windup times."""

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
entries = data["entries"]

N = np.array([e["N_carbon"] for e in entries])
t_spiral = np.array([e["time_spiral_us"] for e in entries])
t_windup = np.array([e["time_windup_us"] for e in entries])

t_spiral_ms = t_spiral / 1000.0
t_windup_ms = t_windup / 1000.0

# --- Plot 1: Spiral computation time ---
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(N, t_spiral_ms, 'o', markersize=4, alpha=0.7)
ax.set_xlabel('N (carbon atoms)')
ax.set_ylabel('Time (ms)')
ax.set_title('Canonical Spiral Computation Time — C28 GC transforms')
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(out_dir, 'spiral-computation-time.png'), dpi=150)
print(f"Saved {os.path.join(out_dir, 'spiral-computation-time.png')}")

# --- Plot 2: Spiral windup time ---
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(N, t_windup_ms, 's', markersize=4, alpha=0.7, color='tab:orange')
ax.set_xlabel('N (carbon atoms)')
ax.set_ylabel('Time (ms)')
ax.set_title('Spiral Windup Time — C28 GC transforms')
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(out_dir, 'spiral-windup-time.png'), dpi=150)
print(f"Saved {os.path.join(out_dir, 'spiral-windup-time.png')}")

# --- Plot 3: Spiral computation time per vertex ---
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(N, t_spiral / N, 'o', markersize=4, alpha=0.7)
ax.set_xlabel('N (carbon atoms)')
ax.set_ylabel('Time per vertex (μs)')
ax.set_title('Canonical Spiral Computation Time per Vertex — C28 GC transforms')
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(out_dir, 'spiral-computation-per-vertex.png'), dpi=150)
print(f"Saved {os.path.join(out_dir, 'spiral-computation-per-vertex.png')}")

# --- Plot 4: Spiral windup time per vertex ---
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(N, t_windup / N, 's', markersize=4, alpha=0.7, color='tab:orange')
ax.set_xlabel('N (carbon atoms)')
ax.set_ylabel('Time per vertex (μs)')
ax.set_title('Spiral Windup Time per Vertex — C28 GC transforms')
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(out_dir, 'spiral-windup-per-vertex.png'), dpi=150)
print(f"Saved {os.path.join(out_dir, 'spiral-windup-per-vertex.png')}")
