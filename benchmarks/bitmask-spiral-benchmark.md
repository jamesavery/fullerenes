# Bitmask spiral benchmark: `remaining_graph` copy vs per-node bitmask

Date: 2026-02-24
Branch: `cleanup_gc`
Compiler: g++ -O2 -std=c++20
Hardware: same machine, back-to-back runs

## Change

Replaced `PlanarGraph remaining_graph(Graph(neighbours,true))` in
`get_spiral_implementation` with a flat `vector<uint16_t> active_mask(N)`.

The old code copied the entire graph (N heap-allocated vectors) on each of
the ~120 calls per canonical spiral search. The new code allocates one flat
array (2N bytes) and uses bit operations to track which neighbors are still
active, referencing the original `this->neighbours` read-only.

## Test case

GC(k,0) transforms of C28 (N_dual = 16k^2), ranging from k=1 (N=16) to
k=56 (N=43,906, corresponding to ~C88,000).

## Results

```
   k   N_dual   Old (ms)   New (ms)   Speedup   Old us/node   New us/node
   1       16      0.10       0.06      1.6x         6.0           3.7
   2       58      0.50       0.19      2.6x         8.5           3.3
   3      128      1.13       0.43      2.6x         8.8           3.4
   4      226      2.05       0.76      2.7x         9.1           3.4
   5      352      3.17       1.20      2.6x         9.0           3.4
   6      506      4.42       1.67      2.6x         8.7           3.3
   8      898      7.79       2.96      2.6x         8.7           3.3
  10    1,402     11.91       4.61      2.6x         8.5           3.3
  15    3,152     26.51       9.95      2.7x         8.4           3.2
  20    5,602     47.64      17.60      2.7x         8.5           3.1
  25    8,752     76.70      30.58      2.5x         8.8           3.5
  30   12,602    110.70      39.44      2.8x         8.8           3.1
  35   17,152    152.24      54.46      2.8x         8.9           3.2
  40   22,402    198.77      71.04      2.8x         8.9           3.2
  45   28,352    251.86      90.18      2.8x         8.9           3.2
  50   35,002    311.38     111.73      2.8x         8.9           3.2
  56   43,906    392.27     138.85      2.8x         8.9           3.2
```

## Summary

- Consistent **2.6-2.8x speedup** across all sizes.
- Both old and new scale linearly (constant us/node), confirming the earlier
  `list -> Deque` change already fixed the O(N^1.5) allocation scaling.
- The bitmask eliminates ~120 x N vector allocations per canonical search,
  reducing us/node from ~8.9 to ~3.2.
