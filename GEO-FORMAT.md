# GEO-FORMAT.md — the compact binary geometry format (`.geo`), version 0

Status: implemented (`include/fullerenes/geo-format.hh`,
`src/c++/geo-format.cc`; tests in `tests/geo-io-test.cc`). The C++ interface
is in §13; the measurements behind the format's choices are in §14.

## 1. Purpose

A `.geo` file is an archive of polyhedra (or of bare point sets, or of bare
surface graphs) with fixed-size records, so that record *i* is one seek away.
It is designed for storing geometry of whole isomer spaces — e.g. the
Alexandrov embeddings of all C<sub>N</sub> isomers — with as much precision per
byte as possible, while staying a general-purpose format:

- **coordinates**: fixed point of any width 2–30 bits, `f32`, `f64`, or none;
- **connectivity** (optional): an oriented cell decomposition of the sphere,
  stored as a half-edge twin matching, so non-simplicial Δ-complexes
  (self-loops, parallel edges) are representable;
- **capacities**: a file-wide vertex capacity and degree range make every
  record the same size, including graphs with varying degrees and records
  with a varying vertex count;
- **integrity**: a 64-bit checksum that is updated in O(record) per append.

## 2. Conventions

- All multi-byte fields are **little-endian**; floats are IEEE-754.
- **Bit streams** are least-significant-bit first: stream bit *k* is bit
  (*k* mod 8) of byte ⌊*k*/8⌋, and a *b*-bit field written at stream position
  *p* has its least significant bit at *p*.
- `bit_length(x)` is the number of binary digits of *x* (`bit_length(0) = 0`);
  a value in [0, *m*) is stored in `bit_length(m−1)` bits.
- Every bit that carries no data — reserved fields, unused slots, padding —
  **is 0**, and readers check it.
- *N* is the header's vertex capacity; *n* ≤ *N* is a record's vertex count.

## 3. File layout

```
[header: 32 bytes] [record 0] [record 1] ... [record count−1] [ignored tail]
```

Every record has the same size *R* (§5.2), derived from the header alone, so
record *i* starts at byte `32 + i·R`. Bytes after record `count−1` are the
remains of an interrupted append (§10); readers ignore them.

## 4. Header

```
off  size  field
0    1     spec       bits 7–5 version (= 0)
                      bit  4   graph: records carry connectivity
                      bit  3   reserved (0)
                      bit  2   offset: records carry an offset (fixed point only)
                      bits 1–0 type: 0 fixed point, 1 no coordinates, 2 f32, 3 f64
1    3     N          vertex capacity, u24, ≥ 1 (≥ 3 with graph)
4    1     width      fixed point: 2..30; otherwise 0
5    1     flags      bit 0  scale per record (fixed point only)
                      bit 1  vertex count per record
                      bit 2  exact triangulation (graph only)
                      bits 3–7 reserved (0)
6    1     deg_min    graph: smallest admissible degree; otherwise 0
7    1     deg_bits   graph: 0..27; otherwise 0
8    4     scale      f32: the file's scale (fixed point without flag bit 0),
                      finite and > 0; otherwise 0
12   4     reserved   0
16   8     count      u64, number of records
24   8     checksum   u64 (§9)
```

Bytes 0–15 never change after the file is created; bytes 16–31 change with
every append. Type 1 (no coordinates) requires the graph bit.

## 5. Records

### 5.1 Layout

Fields in this order, each present only under its condition:

```
[scale     f32]                  type 0 with flag bit 0
[offset    3 × f32]              offset bit
[coords    3N × f32 | f64]       type 2 | 3; x, y, z per vertex slot
bit stream:
  [n         bit_length(N)]      flag bit 1
  [coords    3N × width]         type 0; x, y, z per vertex slot (§6.1)
  [degrees   N × deg_bits]       graph (§7.3)
  [twins     E_cap × tw]         graph (§7.5)
[padding]                        zeros, to a multiple of a
```

with the graph capacities of §7.2 (`E_cap`, `tw = bit_length(2·E_cap − 1)`)
and the record alignment `a` = 8 for type 3, 4 if the record holds any `f32`
(type 2, a per-record scale, or an offset), and 1 otherwise. *R* is a
multiple of *a* and the header is 32 bytes, so every record starts at a
multiple of *a* and every field sits at a multiple of its own size: a
memory-mapped file can be read in place.

### 5.2 Record size

```
P    = 4·[flag bit 0] + 12·[offset bit] + 3N·sizeof(T)·[type 2 or 3]
bits = bit_length(N)·[flag bit 1] + 3N·width·[type 0]
       + (N·deg_bits + E_cap·tw)·[graph bit]
R    = a · ⌈(P + ⌈bits/8⌉) / a⌉
```

## 6. Coordinates

Vertex slots `0..n−1` hold the live vertices; slots `n..N−1` (§8) are zero.

### 6.1 Fixed point (type 0)

With width *w*, bias *B* = 2<sup>w−1</sup> and *q* = 2<sup>w−1</sup> − 1, a
stored *w*-bit value *u* decodes as

    x = o + σ·(u − B)

where σ is the record's scale (flag bit 0) or the header's, and *o* is the
record's offset (offset bit) or (0, 0, 0) — without an offset, coordinates are
assumed to be centred near the origin. Every *u* in [0, 2<sup>w</sup>−1] is a
legal value.

**Writing.** `u = B + llround((x − o)/σ)` (nearest integer, ties away from
zero); a writer only produces |u − B| ≤ q, i.e. u ∈ [1, 2<sup>w</sup>−1], and
refuses a coordinate outside that range (`ValueOutOfRange`, only possible with
a file scale) or a non-finite one.

**Per-record scale.** With *s* = max |x − o| over the record's live
coordinates, σ is the smallest `f32` ≥ *s*/*q*, or 1 when *s* = 0. Then
|x − o|/σ ≤ q, so every coordinate is in range.

**File scale.** Fixed when the file is created, so it must cover every record
ever appended. A proven bound: every edge of length ≤ ℓ<sub>max</sub> gives
|x<sub>u</sub> − x<sub>v</sub>| ≤ d<sub>G</sub>(u, v)·ℓ<sub>max</sub>, so a
record centred on its bounding box has half-extent ≤ diam(G)·ℓ<sub>max</sub>/2,
and diam(G) ≤ N/5 + 1 for a fullerene graph on N vertices (Andova, Došlić,
Krnc, Lužar, Škrekovski, MATCH 68 (2012) 109–130). For isomer spaces that
bound costs 1.6–4× in precision against a per-record scale (§14).

**Offset.** The writer's choice; the bounding-box midpoint rounded to `f32`
minimises *s*. Quantisation uses the stored `f32` value.

**Precision.** |x − x̂| ≤ σ/2 per coordinate, up to one double-precision
rounding in (x − o)/σ and one in the final sum. The product σ·(u − B) is exact
in a double: σ has at most 24 significant bits and |u − B| ≤ 2<sup>29</sup>
because *w* ≤ 30.

### 6.2 Floating point (types 2, 3)

3N little-endian `f32` or `f64` values, byte-aligned at the start of the
record (there is no scale or offset). Values must be finite. `f64` is lossless.

### 6.3 No coordinates (type 1)

The record holds only the bit stream: an archive of surface graphs.

## 7. Graph

### 7.1 Admitted complexes

A record's graph is a connected, closed, oriented combinatorial surface of
genus 0 — a cell decomposition of the sphere — in which every face has at
least 3 sides (counted with multiplicity). Self-loops and parallel edges are
allowed, so non-simplicial intrinsic Delaunay triangulations are
representable. The rotation system is counter-clockwise seen from outside
(the library's orientation invariant), and faces are traced counter-clockwise.

### 7.2 Capacities

Let *A* = Σ deg = 2*E* be the number of half-edges. Euler's formula
n − E + F = 2 together with 2E = Σ (face sides) ≥ 3F gives E ≤ 3n − 6, with
equality iff every face is a triangle; and E ≤ ⌊n·deg<sub>max</sub>/2⌋. Hence
every admitted record fits

    E_cap = min(⌊N·deg_max/2⌋, 3N − 6),   deg_max = deg_min + 2^deg_bits − 1,

and the twin block has `E_cap` slots of `tw = bit_length(2·E_cap − 1)` bits.
A record with *E* < `E_cap` leaves the last `E_cap − E` slots zero.

### 7.3 Degrees

`N` slots of `deg_bits` bits; slot *v* < *n* holds deg(v) − deg<sub>min</sub>.
With `deg_bits = 0` the graph is regular of degree deg<sub>min</sub>. A degree
never exceeds *A* ≤ 6N − 12 < 2<sup>27</sup>, hence `deg_bits` ≤ 27.

### 7.4 Half-edge order

Half-edges are numbered 0..A−1 grouped by origin: first vertex 0's outgoing
half-edges, then vertex 1's, and so on; within a vertex in counter-clockwise
order seen from outside, starting anywhere (the writer uses its in-memory row
start or `v_out`). With off(v) = Σ<sub>u<v</sub> deg(u), vertex *v* owns
half-edges off(v) .. off(v)+deg(v)−1. `rot(h)` is the next half-edge of the
same origin counter-clockwise (wrapping), `rot⁻¹(h)` the previous one.

### 7.5 Twin block

The twin involution is stored once per edge. Walk h = 0, 1, …, A−1; whenever
*h* is still unmatched, the next entry *t* is its twin: it must satisfy
h < t < A with *t* unmatched, and both become matched. This reads exactly
E = A/2 entries, and entry *k* defines **edge k** with half-edges
h<sub>k</sub> < t<sub>k</sub>.

### 7.6 Decoding

- origin(h) is the vertex owning slot *h*; target(h) = origin(twin(h)).
- The face successor is `next(h) = rot⁻¹(twin(h))`: the half-edge before
  twin(h) in target(h)'s rotation. Faces are the cycles of `next`.
- **DelaunayTriangulation**: half-edges 2k and 2k+1 are h<sub>k</sub> and
  t<sub>k</sub> (the library's `twin(h) = h^1`); `he_origin` and `he_next` follow
  from the rules above, `he_face`/`f_he` from the `next` cycles, and
  `v_out[v]` is half-edge off(v). This is the library identity
  `cw(h) = he_next[twin(h)]`. Only connectivity is stored: the metric
  (`he_length`, angles, `v_orig_degree`) is supplied by the caller.
- **Graph / Polyhedron**: row *v* lists target(h) for the half-edges of *v* in
  order, and the RSR twin position of slot (v, i) is twin(h) − off(target(h)).
  Such a record must be simple (no self-loops or parallel edges).

The encoding is about ⅔ of the size of neighbour rows for simple graphs, and
approaches ½ for large *N*, because rows store each edge twice (§14).

### 7.7 Exact triangulation (flag bit 2)

The writer asserts that every face is a triangle, i.e. E = 3n − 6, and readers
check it. Required when reading a record as a `DelaunayTriangulation`.

## 8. Vertex count per record (flag bit 1)

Records then store *n* ∈ [1, N] ([3, N] with a graph) in `bit_length(N)` bits;
without the flag, n = N. Vertex slots n..N−1 — coordinates and degrees — are
zero, and the per-record scale is taken over the live coordinates only.
Nothing else depends on *n* for locating fields, so the record size stays
fixed.

## 9. Checksum

    seed     = XXH3_64bits(header bytes 0–15)
    h_i      = XXH3_64bits_withSeed(le64(i) ‖ record i, seed)     (R + 8 bytes)
    checksum = Σ_{i < count} h_i  mod 2^64                          (0 when empty)

XXH3 is the 64-bit XXH3 of xxHash ≥ 0.8.0, whose output is frozen; a sanity
value is `XXH3_64bits("") = 0x2d06800538d394c2`. The record bytes include their
padding.

- Appending costs one hash of the new record plus an addition.
- The index inside each hash detects swapped records, the seed detects any
  change to the fixed header, and the sum over `count` records detects
  truncation or a wrong count. Random corruption goes undetected with
  probability 2<sup>−64</sup>. The checksum guards against accidents, not
  deliberate tampering.
- Verification is parallel over records; rewriting record *i* in place
  costs one re-hash.
- Merging two files with equal header bytes 0–15 is concatenating their
  record areas and re-hashing the appended records under their new indices.

## 10. Writing and appending

**Creating** a file writes the header with `count = 0` and `checksum = 0`.
A single-polyhedron write derives the capacities tightly from its one
polyhedron, and later appends must then fit them. A file created by appending
must declare what later records have to fit: the degree range (with a graph)
and *N* (with a vertex count per record). Deriving them from the first record
would be a trap: the first C60 isomer's 12-cone iDT has degrees 4–6, the
whole C60 space 3–8.

**Appending** a record:

1. take an exclusive `flock` on the file (one writer per file);
2. read and validate the header; the caller's options must match it
   (`HeaderMismatch`), and the polyhedron must fit it (`CapacityExceeded`,
   `ValueOutOfRange`);
3. write the record at `32 + count·R` and flush (and `fsync` with
   `geo_options::sync`);
4. write the new `count` and `checksum` — bytes 16–31 — in one write, and
   flush (and `fsync`).

Writing a fresh file encodes every record before the stream is touched, so a
refused record leaves an existing file as it was; the records are written
before the header, so an interrupted write leaves a zero header that no reader
accepts.

A crash before step 4 completes leaves `count` and `checksum` describing the
old, consistent state; the partial record is ignored by readers and
overwritten by the next append. Step 4 writes 16 bytes inside the first disk
sector; should it ever be torn, verification reports it.

The `FILE*` must be opened for reading and writing (`"r+b"` or `"w+b"`); a
stream in append mode (`O_APPEND`) is refused, since step 4 would land at the
end of the file.

**`Polyhedron::to_file`** picks the options from the name:

| name | coordinates | graph |
|---|---|---|
| `name.geo` | `f64` | yes |
| `name.f32.geo` | `f32` | yes |
| `name.q<w>.geo` (e.g. `name.q12.geo`) | fixed point, *w* bits, scale per record, no offset | yes |

`from_file` reads record 0 of any `.geo` file.

## 11. Reading and validation

A reader rejects, with `mesh_io_error`:

| condition | code |
|---|---|
| version ≠ 0; a reserved bit or byte ≠ 0; type 1 without graph | `UnsupportedFormat` |
| width, scale, offset, flags or degree fields inconsistent with type and graph bit; N out of range | `UnsupportedFormat` |
| file shorter than `32 + count·R` | `MalformedFile` |
| index ≥ count | `IndexOutOfRange` |
| non-zero padding or unused slot; non-finite float; n outside its range | `MalformedFile` |
| degree sum odd or above 2·E_cap; a twin entry out of range, out of order or already matched | `InvalidTopology` |
| a face with fewer than 3 sides (≠ 3 with flag bit 2); disconnected; n − E + F ≠ 2 | `InvalidTopology` |
| self-loop or parallel edge when reading into `Polyhedron` / `Graph` | `NonSimplicial` |
| a record flagged as not an exact triangulation, read as `DelaunayTriangulation` | `NotATriangulation` |
| a null `FILE*` | `NullFile` |

All record checks are O(R). Reading a record does not verify the whole-file
checksum; `verify_geo` does, and returns false when it differs.

A writer refuses anything a reader would reject — it runs the same record
checks before writing — and in addition:

| condition | code |
|---|---|
| a vertex count, degree or half-edge count outside the file's capacities | `CapacityExceeded` |
| a non-finite coordinate, or one outside the file scale's range or the `f32` range | `ValueOutOfRange` |
| appending with options that differ from the file's header | `HeaderMismatch` |
| a self-loop or repeated neighbour in a `Polyhedron`'s rows | `NonSimplicial` |

## 12. Worked example

A regular tetrahedron with vertices (1,1,1), (1,−1,−1), (−1,1,−1), (−1,−1,1),
written with 8-bit fixed point, a scale per record, no offset, and its graph
flagged as an exact triangulation; one record.

Counter-clockwise rows: 0: [1,2,3], 1: [0,3,2], 2: [0,1,3], 3: [0,2,1]. The
half-edges 0..11 are therefore 0→1, 0→2, 0→3, 1→0, 1→3, 1→2, 2→0, 2→1, 2→3,
3→0, 3→2, 3→1, with twins [3,6,9,0,11,7,1,5,10,2,8,4]. The walk of §7.5 reads
the entries 3, 6, 9, 11, 7, 10. Capacities: deg<sub>max</sub> = 3,
E<sub>cap</sub> = min(6, 6) = 6, tw = bit_length(11) = 4. For *s* = 1, σ is the
smallest `f32` ≥ 1/127, `0x3c010205` ≈ 0.0078740167, so ±1 are stored as
128 ± 127. The bit stream has 12·8 + 6·4 = 120 bits, and
R = 4·⌈(4 + 15)/4⌉ = 20.

```
00  10                        spec: version 0, graph, type 0
01  04 00 00                  N = 4
04  08                        width 8
05  05                        flags: scale per record, exact triangulation
06  03 00                     deg_min 3, deg_bits 0 (regular)
08  00 00 00 00               file scale: none
0c  00 00 00 00               reserved
10  01 00 00 00 00 00 00 00   count 1
18  6e df b3 90 3f 6b 5f cb   checksum 0xcb5f6b3f90b3df6e
20  05 02 01 3c               record 0: σ
24  ff ff ff  ff 01 01        coordinates (1,1,1) (1,−1,−1)
2a  01 ff 01  01 01 ff        coordinates (−1,1,−1) (−1,−1,1)
30  63 b9 a7                  twin entries 3,6 | 9,11 | 7,10
33  00                        padding to 4
```

The header seed is `XXH3_64bits(bytes 00–0f) = 0x19e3510abbbfcdbe`.
Decoding gives ±127·σ = ±1.00000011, within σ/2 of ±1. For half-edge 0→1,
`next` = rot⁻¹(twin) = rot⁻¹(1→0) = 1→2, then rot⁻¹(2→1) = 2→0: the face
(0, 1, 2), whose normal (4, 4, −4) points away from the centroid.

## 13. C++ interface

The codec works on a neutral record — vertex count, coordinates, degrees and
the twin involution in §7.4 order — in `include/fullerenes/geo-format.hh` and
`src/c++/geo-format.cc`, independent of the graph types; `polyhedron-io.cc`
and `delaunay.cc` adapt it. `xxhash.h` (BSD-2) is vendored into `src/contrib`.

```cpp
enum class geo_type : uint8_t { FIXED = 0, NONE = 1, F32 = 2, F64 = 3 };

struct geo_options {
  geo_type type          = geo_type::F64;
  uint8_t  width         = 0;      // FIXED: 2..30
  float    scale         = 0;      // FIXED: 0 = a scale per record, > 0 = the file's scale
  bool     offset        = false;  // FIXED: an offset per record
  bool     graph         = true;
  bool     triangulation = false;  // flag bit 2
  bool     record_n      = false;  // flag bit 1; N is then a capacity
  uint32_t N             = 0;      // vertex capacity; 0 = derive
  int      deg_min       = -1;     // -1 = derive
  int      deg_bits      = -1;     // -1 = derive
  bool     sync          = false;  // fsync after each write step
};

struct geo_header {
  geo_options opt;                 // resolved
  uint64_t count, checksum;
  uint64_t edge_capacity() const;  // E_cap
  uint64_t record_size() const;    // R
  uint64_t record_offset(uint64_t i) const { return 32 + i * record_size(); }
};

struct geo_record {                // one record, independent of the graph types
  int n;
  std::vector<coord3d> x;          // n points, or none (type NONE)
  std::vector<int> degree;         // n degrees, or none (no graph)
  std::vector<int> twin;           // the twin involution in §7.4 half-edge order
};

namespace geo {                    // the codec
  geo_header parse_header(std::span<const uint8_t, 32>);   std::array<uint8_t, 32> header_bytes(const geo_header&);
  geo_header resolve(const geo_options&, std::span<const geo_record>);
  std::vector<uint8_t> encode(const geo_header&, const geo_record&);
  geo_record decode(const geo_header&, std::span<const uint8_t>);
  void check_surface(const geo_record&, bool triangulation);   bool is_simple(const geo_record&);
  geo_header read_header(FILE*);   geo_record read_record(FILE*, uint64_t index, geo_header* = nullptr);
  bool write(FILE*, const geo_options&, std::span<const geo_record>);
  bool append(FILE*, const geo_options&, const geo_record&);
  bool verify(FILE*);
}

// Polyhedron
static geo_header read_geo_header(FILE *file);
static bool       verify_geo(FILE *file);
static Polyhedron from_geo(FILE *file, uint64_t index = 0);
static Polyhedron from_geo(FILE *file, const PlanarGraphView &G, uint64_t index = 0); // files without graph
static bool       to_geo(const Polyhedron &P, FILE *file, bool append = false, const geo_options &opt = {});
static bool       to_geo(std::span<const Polyhedron> Ps, FILE *file, const geo_options &opt = {});

// DelaunayTriangulation: connectivity only, always flagged as an exact triangulation.
// The reader asks the caller for the metric once the topology is built: every edge's
// length (by its even half-edge 2k, edge k of the file) and every vertex's original degree.
static bool to_geo(const DelaunayTriangulation &D, std::span<const coord3d> x, FILE *file,
                   bool append = false, geo_options opt = {});
using GeoEdgeLength = std::function<double(const DelaunayTriangulation &D, int h)>;
using GeoOrigDegree = std::function<int(int v)>;
static DelaunayTriangulation from_geo(FILE *file, uint64_t index, const GeoEdgeLength &length,
                                      const GeoOrigDegree &orig_degree, std::vector<coord3d> *x = nullptr);
```

Python (`src/pybind`, by path): `Polyhedron.to_geo(path, append=False,
options=GeoOptions())`, `Polyhedron.to_geo_batch(path, polyhedra, options)`,
`Polyhedron.from_geo(path, index=0, graph=None)` (`graph`: any bound graph or
geometry object, for files without connectivity), `Polyhedron.read_geo_header`
and `Polyhedron.verify_geo`, with `GeoType`, `GeoOptions` (keyword-only) and
`GeoHeader`. `DelaunayTriangulation.to_geo(path, points=None, append=False,
options)` writes the connectivity; reading takes two steps, since the file has
no metric: `DelaunayTriangulation.read_geo_topology(path, index)` returns
`(he_origin, he_next, points)` in the numbering `from_geo` builds, and
`DelaunayTriangulation.from_geo(path, index, lengths=..., orig_degree=...)`
takes `lengths[k]` for edge k and an original degree per vertex (or one for
all). An unopenable path raises `OSError`, a record index past the end
`IndexError`, invalid options `ValueError`, and the other `mesh_io_error`s
`RuntimeError`.

The derive fields (`N = 0`, `deg_min = deg_bits = -1`) are taken tightly from
the records a file is created with; `to_geo(Ps, …)` takes them over the whole
batch. Caller mistakes — an offset or scale with a float type, a width out of
range, a missing declaration when appending creates the file, a stream in
append mode or not open for reading and writing — throw
`std::invalid_argument`; file and data problems throw `mesh_io_error` with the
codes of §11. The writers return false when a stdio write fails.

## 14. Design notes

Measured with `benchmarks/geo_measure.cc` and `benchmarks/geo_measure_report.py`
on the Alexandrov embeddings of all C60 and C80 isomers and samples at C100,
C150 and C200 (46,680 isomers, each realised as the 12-cone dual and the
20–60-cone cubic polytope; no failures, no non-simplicial iDT). Errors are
worst-case position errors in edge lengths, with a scale per record.

| | 8 bits | 12 bits | 16 bits |
|---|---|---|---|
| dual, C60 | 68 B, 1.7e-2 | 88 B, 1.1e-3 | 104 B, 6.6e-5 |
| dual, C200 | 68 B, 3.4e-2 | 88 B, 2.1e-3 | 104 B, 1.3e-4 |
| cubic, C60 | 412 B, 2.7e-2 | 504 B, 1.7e-3 | 592 B, 1.1e-4 |
| cubic, C200 | 404 B, 5.8e-2 | 496 B, 3.6e-3 | 584 B, 2.3e-4 |

- **Scale per record**: the same worst-case error as a scale fitted to the
  whole batch, a 30% lower median, and 1.6–4.2× lower error than the proven
  file-scale bound; 4 bytes. It also makes overflow impossible, which matters
  because whole isomer spaces are written by appending.
- **Free width**: each 2 bits buy 4× precision for about 8 B (dual) or 44 B
  (cubic) per record, so 8 and 16 bits alone are too coarse. Per-axis scales
  lose at equal size: C60 dual at 14 bits with one scale is 96 B at 2.7e-4,
  at 12 bits with three scales 96 B at 7.2e-4.
- **Twin matching** instead of neighbour rows: C60 cubic 90 B vs 135 B,
  C1000 cubic 2250 B vs 3750 B, C70 dual 110 B vs 167 B, 12-cone iDT 29 B vs
  66 B — and it represents Δ-complexes.
- **Capacities**: padding the cubic embeddings to 60 cones costs +34% at C60
  but +8% at C150–C200, where the data is; the dual's T-bar wastes at most
  9 B. An offset index for variable-length records was the alternative.
- **Connectivity** is 30–45% of a record (27 B of the dual's 88 B), hence
  optional: coordinates alone suit a list of known graphs.
- **Deferred to a later version**: predictive coordinate coding, i.e. storing
  lattice residuals against positions predicted from already-decoded
  neighbours, which makes the needed width independent of the cage's size and
  shape.
