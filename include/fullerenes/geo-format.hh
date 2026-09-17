#pragma once

// The compact binary geometry format (.geo), version 0.  GEO-FORMAT.md is the
// specification; the section numbers below refer to it.
//
// This is the codec.  It works on a geo_record -- vertex count, coordinates, degrees
// and the twin involution in half-edge order (sec. 7.4) -- and knows nothing of the
// graph types; Polyhedron::{to,from}_geo (polyhedron-io.cc) and
// DelaunayTriangulation::{to,from}_geo (delaunay.cc) adapt it.  Data problems throw
// mesh_io_error, caller mistakes std::invalid_argument; the writers return false on
// a failed stdio write.

#include "fullerenes/geometry.hh"
#include "fullerenes/mesh-io-error.hh"

#include <array>
#include <cstdint>
#include <cstdio>
#include <span>
#include <vector>

enum class geo_type : uint8_t { FIXED = 0, NONE = 1, F32 = 2, F64 = 3 };

// What a writer asks for (sec. 13).  N, deg_min and deg_bits at their "derive" values
// are taken tightly from the records a file is created with; a file created by
// appending must declare the degree range (with a graph) and N (with record_n).
struct geo_options {
  geo_type type          = geo_type::F64;
  uint8_t  width         = 0;      // FIXED: 2..30
  float    scale         = 0;      // FIXED: 0 = a scale per record, > 0 = the file's scale
  bool     offset        = false;  // FIXED: an offset per record
  bool     graph         = true;
  bool     triangulation = false;  // every face a triangle; checked on write and read
  bool     record_n      = false;  // a vertex count per record; N is then a capacity
  uint32_t N             = 0;      // vertex capacity; 0 = derive
  int      deg_min       = -1;     // smallest admissible degree; -1 = derive
  int      deg_bits      = -1;     // degree field width; -1 = derive
  bool     sync          = false;  // fsync after each write step (sec. 10)
};

// A file's header (sec. 4): resolved options plus the mutable count and checksum.
struct geo_header {
  geo_options opt;           // resolved: N >= 1; deg_min, deg_bits >= 0 (both 0 without graph)
  uint64_t count    = 0;
  uint64_t checksum = 0;

  static constexpr uint64_t size = 32;

  bool     record_scale() const { return opt.type == geo_type::FIXED && opt.scale == 0; }
  uint64_t deg_max() const;          // deg_min + 2^deg_bits - 1
  uint64_t edge_capacity() const;    // E_cap (sec. 7.2); 0 without graph
  uint64_t twin_bits() const;        // bit_length(2 E_cap - 1)
  uint64_t record_size() const;      // R (sec. 5.2)
  uint64_t record_offset(uint64_t i) const { return size + i * record_size(); }
};

// One record, independent of the graph types.  Half-edges 0..A-1 are grouped by
// origin in vertex order, each vertex's counter-clockwise seen from outside (sec. 7.4);
// twin[h] is the twin of half-edge h.
struct geo_record {
  int                  n = 0;   // vertex count
  std::vector<coord3d> x;       // n points; empty for type NONE
  std::vector<int>     degree;  // n degrees; empty without graph
  std::vector<int>     twin;    // sum(degree) half-edges; empty without graph
};

namespace geo {

// --- Header ---------------------------------------------------------------------
// Serialise a valid header, and parse one with every header check of sec. 11
// (mesh_io_error UnsupportedFormat).
std::array<uint8_t, geo_header::size> header_bytes(const geo_header& H);
geo_header parse_header(std::span<const uint8_t, geo_header::size> bytes);

// Validate `opt` (std::invalid_argument) and resolve its "derive" fields tightly
// from `records` (count and checksum are left 0).
geo_header resolve(const geo_options& opt, std::span<const geo_record> records);

// --- Records --------------------------------------------------------------------
// encode: the R bytes of `r` under `H`, refusing anything decode would reject
//   @throws mesh_io_error CapacityExceeded, ValueOutOfRange, InvalidTopology;
//           std::invalid_argument when r's arrays do not match n
// decode: the record in `bytes` (exactly R of them), with every record check of sec. 11
//   @throws mesh_io_error MalformedFile, InvalidTopology
std::vector<uint8_t> encode(const geo_header& H, const geo_record& r);
geo_record           decode(const geo_header& H, std::span<const uint8_t> bytes);

// --- Graph ----------------------------------------------------------------------
// off(v) of sec. 7.4 for v = 0..n: off[n] = A.
std::vector<int> row_offsets(std::span<const int> degree);
// next(h) = rot^-1(twin(h)), for all half-edges (sec. 7.6).
std::vector<int> face_successors(const geo_record& r);
// Sec. 7.1: a connected genus-0 surface whose faces have >= 3 sides (exactly 3 with
// `triangulation`), given a valid twin involution.  @throws mesh_io_error InvalidTopology
void check_surface(const geo_record& r, bool triangulation);
// No self-loops and no parallel edges.
bool is_simple(const geo_record& r);

// --- Checksum (sec. 9) ----------------------------------------------------------
uint64_t header_seed(const geo_header& H);
uint64_t record_hash(uint64_t seed, uint64_t index, std::span<const uint8_t> record);

// --- Files ----------------------------------------------------------------------
// read_header: the header, and that the file holds all `count` records
// read_record: record `index`; *header receives the header when non-null
//   @throws mesh_io_error NullFile, MalformedFile, IndexOutOfRange, and decode's
geo_header read_header(FILE* file);
geo_record read_record(FILE* file, uint64_t index, geo_header* header = nullptr);

// write:  a fresh file holding `records` (anything already in the stream is replaced)
// append: one record at the end, creating the file when the stream is empty; the
//         stream must be open for reading and writing and not in append mode
//   @post  result == no stdio write failed
//   @throws as encode; mesh_io_error HeaderMismatch (append); std::invalid_argument
bool write(FILE* file, const geo_options& opt, std::span<const geo_record> records);
bool append(FILE* file, const geo_options& opt, const geo_record& r);

// Whether the stored checksum matches the records.
//   @throws as read_header
bool verify(FILE* file);

}  // namespace geo
