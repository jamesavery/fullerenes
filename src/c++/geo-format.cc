// The .geo codec (GEO-FORMAT.md): header, record bit layout, surface checks,
// checksum, and the crash-consistent append protocol.

#include "fullerenes/geo-format.hh"

#include <algorithm>
#include <bit>
#include <cerrno>
#include <cfloat>
#include <climits>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <system_error>

#include <fcntl.h>
#include <sys/file.h>
#include <sys/stat.h>
#include <unistd.h>

#define XXH_INLINE_ALL
#include "../contrib/xxhash/xxhash.h"

using std::string;
using std::to_string;
using Code = mesh_io_error::Code;

namespace {

constexpr uint32_t N_limit = (uint32_t(1) << 24) - 1;

[[noreturn]] void fail(Code code, const string& what) { throw mesh_io_error(code, ".geo: " + what); }
[[noreturn]] void bad_arg(const string& what) { throw std::invalid_argument(".geo: " + what); }

uint64_t bit_length(uint64_t x) { return std::bit_width(x); }

// ---- Little-endian fields ----
void put_le(uint8_t* p, uint64_t v, int nbytes) {
  for (int k = 0; k < nbytes; k++) p[k] = uint8_t(v >> (8 * k));
}
uint64_t get_le(const uint8_t* p, int nbytes) {
  uint64_t v = 0;
  for (int k = 0; k < nbytes; k++) v |= uint64_t(p[k]) << (8 * k);
  return v;
}
void   put_f32(uint8_t* p, float f)   { put_le(p, std::bit_cast<uint32_t>(f), 4); }
void   put_f64(uint8_t* p, double d)  { put_le(p, std::bit_cast<uint64_t>(d), 8); }
float  get_f32(const uint8_t* p)      { return std::bit_cast<float>(uint32_t(get_le(p, 4))); }
double get_f64(const uint8_t* p)      { return std::bit_cast<double>(get_le(p, 8)); }

// ---- Bit streams, least significant bit first (sec. 2); fields of at most 32 bits ----
class bit_writer {
  uint8_t* p;
  const uint8_t* end;
  uint64_t acc = 0;
  unsigned fill = 0;
  void emit() {
    if (p == end) throw std::logic_error(".geo: bit stream overran its record");
    *p++ = uint8_t(acc);
    acc >>= 8;
    fill -= 8;
  }
public:
  bit_writer(uint8_t* p, const uint8_t* end) : p(p), end(end) {}
  void put(uint64_t v, unsigned b) {
    if (b > 32 || (v >> b) != 0) throw std::logic_error(".geo: value does not fit its bit field");
    acc |= v << fill;
    fill += b;
    while (fill >= 8) emit();
  }
  // Pads the last byte with zero bits; returns the first unwritten byte.
  uint8_t* finish() {
    if (fill > 0) { fill += 8 - fill % 8; emit(); }
    return p;
  }
};

class bit_reader {
  const uint8_t* p;
  const uint8_t* end;
  uint64_t acc = 0;
  unsigned fill = 0;
public:
  bit_reader(const uint8_t* p, const uint8_t* end) : p(p), end(end) {}
  uint64_t get(unsigned b) {
    while (fill < b) {
      if (p == end) throw std::logic_error(".geo: bit stream overran its record");
      acc |= uint64_t(*p++) << fill;
      fill += 8;
    }
    const uint64_t v = acc & ((uint64_t(1) << b) - 1);
    acc >>= b;
    fill -= b;
    return v;
  }
  // Checks that the unread bits of the last byte are zero; returns the first unread byte.
  const uint8_t* finish() {
    if (acc != 0) fail(Code::MalformedFile, "non-zero padding bits after the bit stream");
    return p;
  }
};

// ---- Options ----
void validate(const geo_options& o) {
  switch (o.type) {
    case geo_type::FIXED: case geo_type::NONE: case geo_type::F32: case geo_type::F64: break;
    default: bad_arg("unknown coordinate type");
  }
  if (o.type == geo_type::FIXED) {
    if (o.width < 2 || o.width > 30)
      bad_arg("fixed-point width " + to_string(int(o.width)) + " outside 2..30");
    if (!(std::isfinite(o.scale) && o.scale >= 0)) bad_arg("the file scale must be finite and >= 0");
  } else {
    if (o.width != 0) bad_arg("a width needs fixed point");
    if (o.scale != 0) bad_arg("a scale needs fixed point");
    if (o.offset)     bad_arg("an offset needs fixed point");
  }
  if (o.type == geo_type::NONE && !o.graph) bad_arg("a record without coordinates needs a graph");
  if (o.triangulation && !o.graph) bad_arg("the triangulation flag needs a graph");
  if (o.N > N_limit) bad_arg("N above 2^24 - 1");
  if ((o.deg_min < 0) != (o.deg_bits < 0)) bad_arg("declare deg_min and deg_bits together");
  if (o.deg_min >= 0) {
    if (!o.graph && (o.deg_min != 0 || o.deg_bits != 0)) bad_arg("a degree range needs a graph");
    if (o.deg_min > 255) bad_arg("deg_min above 255");
    if (o.deg_bits > 27) bad_arg("deg_bits above 27");
  }
}

// The capacity conditions a header must meet, or an empty string.
string capacity_problem(const geo_header& H) {
  const auto& o = H.opt;
  if (o.N < 1 || o.N > N_limit) return "N = " + to_string(o.N) + " outside 1..2^24-1";
  if (o.graph) {
    if (o.N < 3) return "a graph needs N >= 3";
    if (H.edge_capacity() < 3) return "the degree range admits no surface (edge capacity below 3)";
  }
  return {};
}

void match(const geo_header& H, const geo_options& o) {
  const auto mismatch = [](const string& field) {
    fail(Code::HeaderMismatch, "appending with a " + field + " different from the file's");
  };
  const auto& f = H.opt;
  if (o.type != f.type)                   mismatch("coordinate type");
  if (o.width != f.width)                 mismatch("width");
  if (o.scale != f.scale)                 mismatch("scale");
  if (o.offset != f.offset)               mismatch("offset flag");
  if (o.graph != f.graph)                 mismatch("graph flag");
  if (o.triangulation != f.triangulation) mismatch("triangulation flag");
  if (o.record_n != f.record_n)           mismatch("record_n flag");
  if (o.N != 0 && o.N != f.N)             mismatch("vertex capacity N");
  if (o.deg_min >= 0 && (o.deg_min != f.deg_min || o.deg_bits != f.deg_bits))
    mismatch("degree range");
}

// ---- Half-edges (sec. 7.4) ----
struct half_edges {
  std::vector<int> off;     // off[v], v = 0..n
  std::vector<int> origin;  // origin[h], h = 0..A-1
  explicit half_edges(const geo_record& r) : off(geo::row_offsets(r.degree)), origin(off.back()) {
    for (int v = 0; v < r.n; v++)
      std::fill(origin.begin() + off[v], origin.begin() + off[v + 1], v);
  }
  int A() const { return off.back(); }
  // The half-edge before h in its origin's counter-clockwise rotation.
  int rot_prev(int h) const {
    const int v = origin[h];
    return h == off[v] ? off[v + 1] - 1 : h - 1;
  }
};

// A fixed-point-free involution on the record's half-edges.
void check_twins(const geo_record& r, int A) {
  if (r.twin.size() != size_t(A))
    bad_arg("twin has " + to_string(r.twin.size()) + " entries for " + to_string(A) + " half-edges");
  for (int h = 0; h < A; h++) {
    const int t = r.twin[h];
    if (t < 0 || t >= A || t == h || r.twin[t] != h)
      fail(Code::InvalidTopology, "twin is not a fixed-point-free involution at half-edge " + to_string(h));
  }
}

// ---- Files ----
FILE* checked(FILE* file) {
  if (!file) fail(Code::NullFile, "null FILE*");
  return file;
}

uint64_t file_size(FILE* file) {
  struct stat st;
  if (fstat(fileno(file), &st) != 0) fail(Code::MalformedFile, "cannot stat the stream");
  return uint64_t(st.st_size);
}

void read_exact(FILE* file, uint64_t offset, uint8_t* buf, size_t n, const char* what) {
  if (offset > uint64_t(std::numeric_limits<off_t>::max()) || fseeko(file, off_t(offset), SEEK_SET) != 0
      || fread(buf, 1, n, file) != n)
    fail(Code::MalformedFile, string("truncated file: cannot read ") + what);
}

bool write_exact(FILE* file, uint64_t offset, const uint8_t* buf, size_t n) {
  return offset <= uint64_t(std::numeric_limits<off_t>::max())
      && fseeko(file, off_t(offset), SEEK_SET) == 0 && fwrite(buf, 1, n, file) == n;
}

bool flush(FILE* file, bool sync) {
  return fflush(file) == 0 && (!sync || fsync(fileno(file)) == 0);
}

void require_writable(FILE* file, bool need_read) {
  const int flags = fcntl(fileno(checked(file)), F_GETFL);
  if (flags < 0) bad_arg("cannot query the stream's mode");
  if (flags & O_APPEND) bad_arg("the stream is in append mode; open it with \"r+b\" or \"w+b\"");
  const int mode = flags & O_ACCMODE;
  if (need_read ? mode != O_RDWR : mode == O_RDONLY)
    bad_arg(need_read ? "appending needs a stream open for reading and writing"
                      : "the stream is not open for writing");
}

// One writer per file (sec. 10).
class exclusive_lock {
  int fd;
public:
  explicit exclusive_lock(FILE* file) : fd(fileno(file)) {
    if (::flock(fd, LOCK_EX) != 0) throw std::system_error(errno, std::generic_category(), ".geo: flock");
  }
  ~exclusive_lock() { ::flock(fd, LOCK_UN); }
  exclusive_lock(const exclusive_lock&) = delete;
  exclusive_lock& operator=(const exclusive_lock&) = delete;
};

}  // namespace

// ================================================================================
// Header
// ================================================================================

uint64_t geo_header::deg_max() const {
  return uint64_t(opt.deg_min) + (uint64_t(1) << opt.deg_bits) - 1;
}

uint64_t geo_header::edge_capacity() const {
  if (!opt.graph || opt.N < 3) return 0;
  const uint64_t N = opt.N;
  return std::min(N * deg_max() / 2, 3 * N - 6);
}

uint64_t geo_header::twin_bits() const {
  const uint64_t E = edge_capacity();
  return E > 0 ? bit_length(2 * E - 1) : 0;
}

uint64_t geo_header::record_size() const {
  const uint64_t N = opt.N;
  const uint64_t prefix = (record_scale() ? 4 : 0) + (opt.offset ? 12 : 0)
                        + (opt.type == geo_type::F32 ? 12 * N : opt.type == geo_type::F64 ? 24 * N : 0);
  const uint64_t bits = (opt.record_n ? bit_length(N) : 0)
                      + (opt.type == geo_type::FIXED ? 3 * N * opt.width : 0)
                      + (opt.graph ? N * opt.deg_bits + edge_capacity() * twin_bits() : 0);
  const uint64_t align = opt.type == geo_type::F64 ? 8 : prefix > 0 ? 4 : 1;
  const uint64_t bytes = prefix + (bits + 7) / 8;
  return (bytes + align - 1) / align * align;
}

std::array<uint8_t, geo_header::size> geo::header_bytes(const geo_header& H) {
  const auto& o = H.opt;
  std::array<uint8_t, geo_header::size> b{};
  b[0] = uint8_t((o.graph ? 0x10 : 0) | (o.offset ? 0x04 : 0) | uint8_t(o.type));
  put_le(&b[1], o.N, 3);
  b[4] = o.width;
  b[5] = uint8_t((H.record_scale() ? 1 : 0) | (o.record_n ? 2 : 0) | (o.triangulation ? 4 : 0));
  b[6] = uint8_t(o.deg_min);
  b[7] = uint8_t(o.deg_bits);
  put_f32(&b[8], o.type == geo_type::FIXED && !H.record_scale() ? o.scale : 0.0f);
  put_le(&b[16], H.count, 8);
  put_le(&b[24], H.checksum, 8);
  return b;
}

geo_header geo::parse_header(std::span<const uint8_t, geo_header::size> b) {
  const auto bad = [](const string& what) { fail(Code::UnsupportedFormat, what); };
  geo_header H;
  auto& o = H.opt;

  const uint8_t spec = b[0];
  if (spec >> 5) bad("version " + to_string(spec >> 5) + "; this reader knows version 0");
  if (spec & 0x08) bad("reserved bit 3 of the spec byte is set");
  o.graph  = spec & 0x10;
  o.offset = spec & 0x04;
  o.type   = geo_type(spec & 0x03);
  o.N      = uint32_t(get_le(&b[1], 3));
  o.width  = b[4];

  const uint8_t flags = b[5];
  if (flags & 0xF8) bad("reserved flag bits are set");
  const bool record_scale = flags & 1;
  o.record_n      = flags & 2;
  o.triangulation = flags & 4;
  o.deg_min  = b[6];
  o.deg_bits = b[7];
  const uint32_t scale_bits = uint32_t(get_le(&b[8], 4));
  if (get_le(&b[12], 4) != 0) bad("reserved header bytes 12-15 are not zero");

  const bool fixed = o.type == geo_type::FIXED;
  if (o.type == geo_type::NONE && !o.graph) bad("type 1 (no coordinates) without a graph");
  if (fixed ? (o.width < 2 || o.width > 30) : o.width != 0)
    bad("width " + to_string(int(o.width)) + " for coordinate type " + to_string(int(o.type)));
  if (!fixed && (record_scale || o.offset)) bad("a scale or offset flag without fixed point");
  if (o.triangulation && !o.graph) bad("the triangulation flag without a graph");
  if (!o.graph && (o.deg_min != 0 || o.deg_bits != 0)) bad("degree fields without a graph");
  if (o.deg_bits > 27) bad("deg_bits " + to_string(o.deg_bits) + " above 27");
  if (fixed && !record_scale) {
    o.scale = std::bit_cast<float>(scale_bits);
    if (!(std::isfinite(o.scale) && o.scale > 0)) bad("the file scale is not finite and positive");
  } else if (scale_bits != 0) {
    bad("a file scale in a file without one");
  }
  if (const string p = capacity_problem(H); !p.empty()) bad(p);

  H.count    = get_le(&b[16], 8);
  H.checksum = get_le(&b[24], 8);
  return H;
}

geo_header geo::resolve(const geo_options& opt, std::span<const geo_record> records) {
  validate(opt);
  geo_header H;
  H.opt = opt;
  auto& o = H.opt;

  if (o.N == 0) {
    if (records.empty()) bad_arg("N must be declared for a file without records");
    int n = 0;
    for (const auto& r : records) n = std::max(n, r.n);
    if (n < 1 || uint64_t(n) > N_limit)
      fail(Code::CapacityExceeded, "vertex count " + to_string(n) + " outside 1..2^24-1");
    o.N = uint32_t(n);
  }

  if (!o.graph) {
    o.deg_min = o.deg_bits = 0;
  } else if (o.deg_min < 0) {
    if (records.empty()) bad_arg("the degree range must be declared for a file without records");
    int lo = INT_MAX, hi = INT_MIN;
    for (const auto& r : records)
      for (int d : r.degree) { lo = std::min(lo, d); hi = std::max(hi, d); }
    if (lo > hi) bad_arg("the records carry no degrees");
    if (lo < 0 || lo > 255)
      fail(Code::CapacityExceeded, "smallest degree " + to_string(lo) + " outside 0..255");
    o.deg_min  = lo;
    o.deg_bits = int(bit_length(uint64_t(hi - lo)));
    if (o.deg_bits > 27) fail(Code::CapacityExceeded, "degree range " + to_string(lo) + ".." + to_string(hi));
  }

  if (const string p = capacity_problem(H); !p.empty()) bad_arg(p);
  return H;
}

// ================================================================================
// Graph
// ================================================================================

std::vector<int> geo::row_offsets(std::span<const int> degree) {
  std::vector<int> off(degree.size() + 1, 0);
  for (size_t v = 0; v < degree.size(); v++) off[v + 1] = off[v] + degree[v];
  return off;
}

std::vector<int> geo::face_successors(const geo_record& r) {
  const half_edges he(r);
  std::vector<int> next(he.A());
  for (int h = 0; h < he.A(); h++) next[h] = he.rot_prev(r.twin[h]);
  return next;
}

void geo::check_surface(const geo_record& r, bool triangulation) {
  const half_edges he(r);
  const int A = he.A();
  for (int v = 0; v < r.n; v++)
    if (r.degree[v] == 0) fail(Code::InvalidTopology, "vertex " + to_string(v) + " has no edge");

  // Faces: the cycles of next = rot^-1 . twin, a permutation.
  std::vector<char> seen(A, 0);
  int64_t F = 0;
  for (int h = 0; h < A; h++) {
    if (seen[h]) continue;
    int sides = 0, g = h;
    do {
      seen[g] = 1;
      g = he.rot_prev(r.twin[g]);
      if (++sides > A) throw std::logic_error(".geo: a face cycle did not close");
    } while (g != h);
    if (sides < 3)
      fail(Code::InvalidTopology, "a face with " + to_string(sides) + " sides (at least 3 required)");
    if (triangulation && sides != 3)
      fail(Code::InvalidTopology, "a " + to_string(sides) + "-sided face in a record flagged as a triangulation");
    F++;
  }

  // Connectivity, from vertex 0 along the half-edges.
  std::vector<char> reached(r.n, 0);
  std::vector<int> queue{0};
  reached[0] = 1;
  for (size_t i = 0; i < queue.size(); i++)
    for (int h = he.off[queue[i]]; h < he.off[queue[i] + 1]; h++) {
      const int w = he.origin[r.twin[h]];
      if (!reached[w]) { reached[w] = 1; queue.push_back(w); }
    }
  if (queue.size() != size_t(r.n))
    fail(Code::InvalidTopology, "the graph is disconnected");

  const int64_t chi = int64_t(r.n) - A / 2 + F;
  if (chi != 2)
    fail(Code::InvalidTopology, "Euler characteristic " + to_string(chi) + " (a sphere has 2)");
}

bool geo::is_simple(const geo_record& r) {
  const half_edges he(r);
  std::vector<int> last(r.n, -1);   // last[w] == v: v already has an edge to w
  for (int v = 0; v < r.n; v++)
    for (int h = he.off[v]; h < he.off[v + 1]; h++) {
      const int w = he.origin[r.twin[h]];
      if (w == v || last[w] == v) return false;
      last[w] = v;
    }
  return true;
}

// ================================================================================
// Records
// ================================================================================

std::vector<uint8_t> geo::encode(const geo_header& H, const geo_record& r) {
  const auto& o = H.opt;
  const uint64_t N = o.N;
  const bool fixed = o.type == geo_type::FIXED;

  const int n_min = o.graph ? 3 : 1;
  if (r.n < n_min || uint64_t(r.n) > N)
    fail(Code::CapacityExceeded, "vertex count " + to_string(r.n) + " outside "
                                 + to_string(n_min) + ".." + to_string(N));
  if (!o.record_n && uint64_t(r.n) != N)
    fail(Code::CapacityExceeded, to_string(r.n) + " vertices in a file whose records have exactly "
                                 + to_string(N) + " (a capacity needs record_n)");
  const size_t n_coords = o.type == geo_type::NONE ? 0 : size_t(r.n);
  if (r.x.size() != n_coords)
    bad_arg(to_string(r.x.size()) + " coordinates for " + to_string(n_coords));
  for (const auto& p : r.x)
    for (int k = 0; k < 3; k++)
      if (!std::isfinite(p[k])) fail(Code::ValueOutOfRange, "a non-finite coordinate");

  int A = 0;
  if (o.graph) {
    if (r.degree.size() != size_t(r.n))
      bad_arg(to_string(r.degree.size()) + " degrees for " + to_string(r.n) + " vertices");
    uint64_t sum = 0;
    for (int v = 0; v < r.n; v++) {
      const int d = r.degree[v];
      if (d < o.deg_min || uint64_t(d) > H.deg_max())
        fail(Code::CapacityExceeded, "vertex " + to_string(v) + " has degree " + to_string(d)
                                     + " outside the file's " + to_string(o.deg_min) + ".." + to_string(H.deg_max()));
      sum += uint64_t(d);
    }
    if (sum > 2 * H.edge_capacity())
      fail(Code::CapacityExceeded, to_string(sum) + " half-edges exceed the file's capacity of "
                                   + to_string(2 * H.edge_capacity()));
    A = int(sum);
    check_twins(r, A);
    check_surface(r, o.triangulation);
  } else if (!r.degree.empty() || !r.twin.empty()) {
    bad_arg("graph arrays for a file without graph");
  }

  std::vector<uint8_t> out(H.record_size(), 0);
  uint8_t* p = out.data();
  const uint8_t* const end = out.data() + out.size();

  // Fixed-point frame (sec. 6.1).
  double origin[3] = {0, 0, 0};
  double sigma = o.scale;
  if (fixed) {
    float offset[3] = {0, 0, 0};
    if (o.offset)
      for (int k = 0; k < 3; k++) {
        double lo = r.x[0][k], hi = r.x[0][k];
        for (const auto& q : r.x) { lo = std::min(lo, q[k]); hi = std::max(hi, q[k]); }
        offset[k] = float(lo / 2 + hi / 2);
        if (!std::isfinite(offset[k])) fail(Code::ValueOutOfRange, "the offset exceeds the f32 range");
        origin[k] = offset[k];
      }
    if (H.record_scale()) {
      const double q = double((uint64_t(1) << (o.width - 1)) - 1);
      double s = 0;
      for (const auto& x : r.x)
        for (int k = 0; k < 3; k++) s = std::max(s, std::fabs(x[k] - origin[k]));
      float f = 1.0f;
      if (s > 0) {
        const double t = s / q;
        f = float(t);
        if (double(f) < t) f = std::nextafter(f, float(INFINITY));
        if (!std::isfinite(f)) fail(Code::ValueOutOfRange, "the coordinates exceed an f32 scale");
      }
      put_f32(p, f);
      p += 4;
      sigma = f;
    }
    if (o.offset)
      for (int k = 0; k < 3; k++) { put_f32(p, offset[k]); p += 4; }
  } else if (o.type == geo_type::F32 || o.type == geo_type::F64) {
    for (uint64_t v = 0; v < N; v++)
      for (int k = 0; k < 3; k++) {
        const double x = v < uint64_t(r.n) ? r.x[v][k] : 0.0;
        if (o.type == geo_type::F32) {
          if (!(std::fabs(x) <= FLT_MAX)) fail(Code::ValueOutOfRange, "a coordinate beyond the f32 range");
          put_f32(p, float(x));
          p += 4;
        } else {
          put_f64(p, x);
          p += 8;
        }
      }
  }

  bit_writer bits(p, end);
  if (o.record_n) bits.put(uint64_t(r.n), unsigned(bit_length(N)));
  if (fixed) {
    const int64_t bias = int64_t(1) << (o.width - 1), q = bias - 1;
    for (uint64_t v = 0; v < N; v++)
      for (int k = 0; k < 3; k++) {
        if (v >= uint64_t(r.n)) { bits.put(0, o.width); continue; }
        const double t = (r.x[v][k] - origin[k]) / sigma;
        const int64_t steps = std::fabs(t) <= double(q) + 0.5 ? std::llround(t) : q + 1;
        if (steps < -q || steps > q)
          fail(Code::ValueOutOfRange, "coordinate " + to_string(r.x[v][k]) + " of vertex " + to_string(v)
                                      + " outside the range of the file scale");
        bits.put(uint64_t(bias + steps), o.width);
      }
  }
  if (o.graph) {
    for (uint64_t v = 0; v < N; v++)
      bits.put(v < uint64_t(r.n) ? uint64_t(r.degree[v] - o.deg_min) : 0, unsigned(o.deg_bits));
    const unsigned tw = unsigned(H.twin_bits());
    std::vector<char> matched(A, 0);
    uint64_t E = 0;
    for (int h = 0; h < A; h++) {
      if (matched[h]) continue;
      const int t = r.twin[h];
      bits.put(uint64_t(t), tw);
      matched[h] = matched[t] = 1;
      E++;
    }
    for (; E < H.edge_capacity(); E++) bits.put(0, tw);
  }
  bits.finish();
  return out;
}

geo_record geo::decode(const geo_header& H, std::span<const uint8_t> bytes) {
  const auto& o = H.opt;
  const uint64_t N = o.N;
  if (bytes.size() != H.record_size()) throw std::logic_error(".geo: decode needs exactly one record");
  const auto malformed = [](const string& what) { fail(Code::MalformedFile, what); };

  const uint8_t* p = bytes.data();
  const uint8_t* const end = p + bytes.size();

  double sigma = o.scale;
  double origin[3] = {0, 0, 0};
  if (H.record_scale()) {
    const float f = get_f32(p);
    p += 4;
    if (!(std::isfinite(f) && f > 0)) malformed("a record scale that is not finite and positive");
    sigma = f;
  }
  if (o.offset)
    for (int k = 0; k < 3; k++) {
      const float f = get_f32(p);
      p += 4;
      if (!std::isfinite(f)) malformed("a non-finite offset");
      origin[k] = f;
    }
  const uint8_t* const floats = p;
  const uint64_t float_size = o.type == geo_type::F32 ? 4 : o.type == geo_type::F64 ? 8 : 0;
  p += 3 * N * float_size;

  bit_reader bits(p, end);
  geo_record r;
  r.n = int(N);
  if (o.record_n) {
    const uint64_t n = bits.get(unsigned(bit_length(N)));
    if (n < (o.graph ? 3u : 1u) || n > N) malformed("vertex count " + to_string(n) + " out of range");
    r.n = int(n);
  }
  const uint64_t n = uint64_t(r.n);

  if (float_size > 0) {
    r.x.resize(n);
    for (uint64_t v = 0; v < N; v++)
      for (int k = 0; k < 3; k++) {
        const uint8_t* q = floats + (3 * v + k) * float_size;
        if (v < n) {
          r.x[v][k] = float_size == 4 ? double(get_f32(q)) : get_f64(q);
          if (!std::isfinite(r.x[v][k])) malformed("a non-finite coordinate");
        } else if (std::any_of(q, q + float_size, [](uint8_t c) { return c != 0; })) {
          malformed("a non-zero unused vertex slot " + to_string(v));
        }
      }
  }
  if (o.type == geo_type::FIXED) {
    r.x.resize(n);
    const int64_t bias = int64_t(1) << (o.width - 1);
    for (uint64_t v = 0; v < N; v++)
      for (int k = 0; k < 3; k++) {
        const uint64_t u = bits.get(o.width);
        if (v < n) r.x[v][k] = origin[k] + sigma * double(int64_t(u) - bias);
        else if (u != 0) malformed("a non-zero unused vertex slot " + to_string(v));
      }
  }

  if (o.graph) {
    r.degree.resize(n);
    uint64_t A = 0;
    for (uint64_t v = 0; v < N; v++) {
      const uint64_t d = bits.get(unsigned(o.deg_bits));
      if (v < n) { r.degree[v] = int(o.deg_min + d); A += uint64_t(r.degree[v]); }
      else if (d != 0) malformed("a non-zero degree slot " + to_string(v) + " beyond the vertex count");
    }
    const uint64_t E_cap = H.edge_capacity();
    if (A % 2 != 0 || A > 2 * E_cap)
      fail(Code::InvalidTopology, "degree sum " + to_string(A) + " is odd or exceeds "
                                  + to_string(2 * E_cap));
    r.twin.assign(A, -1);
    const unsigned tw = unsigned(H.twin_bits());
    uint64_t E = 0;
    for (uint64_t h = 0; h < A; h++) {
      if (r.twin[h] >= 0) continue;
      const uint64_t t = bits.get(tw);
      if (t <= h || t >= A || r.twin[t] >= 0)
        fail(Code::InvalidTopology, "twin entry " + to_string(E) + " = " + to_string(t) + " for half-edge "
                                    + to_string(h) + " (needs h < t < A with t unmatched)");
      r.twin[h] = int(t);
      r.twin[t] = int(h);
      E++;
    }
    for (; E < E_cap; E++)
      if (bits.get(tw) != 0) malformed("a non-zero unused twin slot");
    check_surface(r, o.triangulation);
  }

  for (const uint8_t* q = bits.finish(); q < end; q++)
    if (*q != 0) malformed("non-zero padding bytes");
  return r;
}

// ================================================================================
// Checksum
// ================================================================================

uint64_t geo::header_seed(const geo_header& H) {
  const auto b = header_bytes(H);
  return XXH3_64bits(b.data(), 16);
}

uint64_t geo::record_hash(uint64_t seed, uint64_t index, std::span<const uint8_t> record) {
  XXH3_state_t state;
  XXH3_INITSTATE(&state);
  uint8_t le[8];
  put_le(le, index, 8);
  if (XXH3_64bits_reset_withSeed(&state, seed) != XXH_OK || XXH3_64bits_update(&state, le, 8) != XXH_OK
      || XXH3_64bits_update(&state, record.data(), record.size()) != XXH_OK)
    throw std::logic_error(".geo: XXH3 streaming failed");
  return XXH3_64bits_digest(&state);
}

// ================================================================================
// Files
// ================================================================================

geo_header geo::read_header(FILE* file) {
  checked(file);
  std::array<uint8_t, geo_header::size> b;
  read_exact(file, 0, b.data(), b.size(), "the header");
  const geo_header H = parse_header(b);
  const uint64_t R = H.record_size();
  if (H.count > (std::numeric_limits<uint64_t>::max() - geo_header::size) / R
      || file_size(file) < H.record_offset(H.count))
    fail(Code::MalformedFile, "the file is shorter than its " + to_string(H.count) + " records");
  return H;
}

geo_record geo::read_record(FILE* file, uint64_t index, geo_header* header) {
  const geo_header H = read_header(file);
  if (index >= H.count)
    fail(Code::IndexOutOfRange, "record " + to_string(index) + " of " + to_string(H.count));
  std::vector<uint8_t> bytes(H.record_size());
  read_exact(file, H.record_offset(index), bytes.data(), bytes.size(), "a record");
  if (header) *header = H;
  return decode(H, bytes);
}

bool geo::write(FILE* file, const geo_options& opt, std::span<const geo_record> records) {
  require_writable(file, false);
  geo_header H = resolve(opt, records);
  const uint64_t seed = header_seed(H);

  // Encode everything before touching the file, so a refused record leaves it as it was.
  std::vector<uint8_t> body;
  body.reserve(records.size() * H.record_size());
  for (size_t i = 0; i < records.size(); i++) {
    const auto bytes = encode(H, records[i]);
    H.checksum += record_hash(seed, i, bytes);
    body.insert(body.end(), bytes.begin(), bytes.end());
  }
  H.count = records.size();

  // Records first, header last: a crash midway leaves a zero header, which no reader accepts.
  const auto hb = header_bytes(H);
  return fflush(file) == 0 && ftruncate(fileno(file), 0) == 0
      && write_exact(file, geo_header::size, body.data(), body.size())
      && write_exact(file, 0, hb.data(), hb.size()) && flush(file, opt.sync);
}

bool geo::append(FILE* file, const geo_options& opt, const geo_record& r) {
  require_writable(file, true);
  validate(opt);
  const exclusive_lock lock(file);
  if (fflush(file) != 0) return false;

  geo_header H;
  const bool create = file_size(file) == 0;
  if (create) {
    if (opt.graph && opt.deg_min < 0)
      bad_arg("a file created by appending must declare its degree range");
    if (opt.record_n && opt.N == 0)
      bad_arg("a file created by appending with record_n must declare N");
    H = resolve(opt, std::span<const geo_record>(&r, 1));
  } else {
    H = read_header(file);
    match(H, opt);
  }
  const auto bytes = encode(H, r);

  if (create) {
    const auto hb = header_bytes(H);
    if (!write_exact(file, 0, hb.data(), hb.size()) || !flush(file, opt.sync)) return false;
  }
  // Sec. 10: the record, then count and checksum in one write.
  if (!write_exact(file, H.record_offset(H.count), bytes.data(), bytes.size()) || !flush(file, opt.sync))
    return false;
  H.checksum += record_hash(header_seed(H), H.count, bytes);
  H.count += 1;
  const auto hb = header_bytes(H);
  return write_exact(file, 16, hb.data() + 16, 16) && flush(file, opt.sync);
}

bool geo::verify(FILE* file) {
  const geo_header H = read_header(file);
  const uint64_t seed = header_seed(H);
  std::vector<uint8_t> bytes(H.record_size());
  uint64_t sum = 0;
  for (uint64_t i = 0; i < H.count; i++) {
    if (i == 0) read_exact(file, H.record_offset(0), bytes.data(), bytes.size(), "a record");
    else if (fread(bytes.data(), 1, bytes.size(), file) != bytes.size())
      fail(Code::MalformedFile, "truncated file: cannot read a record");
    sum += record_hash(seed, i, bytes);
  }
  return sum == H.checksum;
}
