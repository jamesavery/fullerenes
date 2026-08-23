#ifndef ISOMERDB_HH
# define ISOMERDB_HH
#include <array>
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <errno.h>

#include "fullerenegraph.hh"
#include "symmetry.hh"

using namespace std;

// The fullerene isomer database ("PDB" = Peter's Database: the text format
// of the legacy Fortran program by P. Schwerdtfeger et al.), one file per
// (N, corpus) holding every isomer with its point group, canonical ring-spiral
// pentagon indices, pentagon/hexagon neighbour indices, Hueckel HOMO data,
// Hamilton-cycle count and NMR pattern.  This class is the C++ reader/writer
// of that format and replaces the Fortran Printdatabase / CompressDatabase
// routines (isomer.f, util.f).
class IsomerDB {
public:
  int N, Nisomers;
  bool IPR, with_ncycham;

  // The three corpora of the database root: All/ (every isomer of C_N),
  // IPR/ (the isolated-pentagon isomers), Nontrivial/ (the isomers with a
  // non-trivial point group).
  enum class Corpus { All, IPR, Nontrivial };

  // Canonical ring-spiral pentagon positions, 1-based, ascending -- the
  // database's primary key.
  using RSPI = std::array<u_int8_t,12>;

  // A field of a record, as a set: which columns a producer fills, and which
  // ones two records differ in.
  //
  // The producer masks name where a column comes from, not what it is:
  // FromDual needs only the dual graph, FromPointGroup the Symmetry class,
  // FromSpectrum a Hueckel analysis, Ncycham a Hamilton-cycle count (minutes
  // per isomer at C110, so it is never implied).
  enum class Field : uint32_t {
    None = 0,
    Group = 1u<<0, RSPI = 1u<<1, PNI = 1u<<2, HNI = 1u<<3,
    NeHOMO = 1u<<4, NedgeHOMO = 1u<<5, HLgap = 1u<<6, Ncycham = 1u<<7, INMR = 1u<<8,
    FromDual        = RSPI | PNI | HNI,
    FromPointGroup  = Group | INMR,
    FromSpectrum    = NeHOMO | NedgeHOMO | HLgap,
    All             = FromDual | FromPointGroup | FromSpectrum | Ncycham,
  };

  struct EntryResult;          // what Entry::from_dual produced, or why not (below)

  struct Entry {
    char group[3];		// point group, right-justified in 3 chars (" Ih", "C2v")
    u_int8_t RSPI[12];		// canonical ring-spiral pentagon positions, 1-based, ascending
    u_int8_t PNI[5];		// PNI[k], k=0..4: #pentagons with k pentagon neighbours (PNI(5) = 12 - sum)
    u_int8_t HNI[6];		// HNI[k], k=0..5: #hexagons with k hexagon neighbours (HNI(6) = N/2-10 - sum)
    u_int8_t NeHOMO;		// electrons in the HOMO level
    u_int8_t NedgeHOMO;		// degeneracy of the HOMO level
    float HLgap;		// HOMO-LUMO gap in units of |beta| (0 for open shells)
    int ncycham;		// number of Hamilton cycles (0 when the file has none: with_ncycham false)
    u_int8_t INMR[6];		// NMR pattern: INMR[2i] orbits of size INMR[2i+1], i=0..2 (unused pairs 0)

    // --- the record as a value ---
    IsomerDB::RSPI rspi() const { IsomerDB::RSPI r; std::copy(RSPI, RSPI+12, r.begin()); return r; }
    PointGroup point_group() const { return PointGroup(string(group, 3)); }
    // The stored gap is the exact gap rounded to the F7.5 column, so two
    // records agree when the column cannot tell them apart: half an F7.5 ulp
    // (5e-6, attained -- two C60 isomers sit exactly on a rounding boundary
    // and round the other way), plus the float half-ulps of the two stored
    // values (<= 1.2e-7 for a gap < 3), plus the eigensolver's own ~1e-8.
    static constexpr float HLgap_column_tol = 5.2e-6f;
    // The largest Hamilton-cycle count the text format's I7 column can hold.
    static constexpr int   ncycham_max = 9999999;

    // The fields of `fields` in which this record and b differ.  Gaps count as
    // equal when they are within HLgap_tol -- exact by default, since two
    // records READ from the database carry the same rounded column; pass
    // HLgap_column_tol when one side was computed.
    Field diff(const Entry& b, Field fields = Field::All, float HLgap_tol = 0.f) const;
    string to_string() const;
    static string to_string_rspi(const IsomerDB::RSPI& r);   // 12 face positions, I4 each

    // The record of the isomer whose dual is D: the fields in `producers` are
    // computed (RSPI always -- it is the key), the others left zero.  The
    // canonical RSPI is FullereneDualView::regular_rspi, the database's key.
    // Definition staged in claude-projects/unfortran until HueckelAnalysis and
    // hamiltonian_cycle_count are promoted.
    // @anchor isomer-record-from-dual
    static EntryResult from_dual(const FullereneDualView& D, Field producers);
    // A filled HOMO level: the Fortran's closed/open distinction, and the
    // condition under which HLgap is a gap rather than the stored 0.
    bool closed_shell() const { return 2*NedgeHOMO == NeHOMO; }
    // An isolated-pentagon isomer: no pentagon has a pentagon neighbour, and
    // (equivalently) no hexagon has fewer than 3 hexagon neighbours.  The IPR
    // file layout stores only the columns this leaves free.
    bool is_ipr() const { return PNI[0] == 12 && !HNI[0] && !HNI[1] && !HNI[2]; }
    // The invariants of a database record of a C_N isomer in a file with
    // these flags: RSPI 12 face positions in [1, N/2+2] strictly ascending;
    // PNI(k), k<5, summing to at most the 12 pentagons, and PNI = (12,0,0,0,0)
    // with no hexagon below 3 hexagon neighbours in an IPR file; HNI(k), k<6,
    // summing to at most the N/2-10 hexagons; NMR orbits accounting for every
    // atom; 1 <= NeHOMO <= 2*NedgeHOMO with NedgeHOMO >= 1; HLgap >= 0 and
    // zero exactly for an open shell (the Fortran convention); a Hamilton
    // count present exactly when the file carries them.  They hold over the
    // entire corpus.  Empty string = valid; otherwise the first violation,
    // for the caller to place in its own context.
    // @anchor isomer-record-invariants
    string invariant_violation(int N, bool IPR, bool with_ncycham) const;
    // The spiral as the graph constructors take it (FullereneGraph(N, rspi),
    // FullereneDual(N, rspi)): 0-based positions.
    vector<int> rspi_zero_based() const { vector<int> r(RSPI, RSPI+12); for(int& x: r) x--; return r; }
    // The neighbour indices as stored: PNI k<5, HNI k<6; the dropped top
    // bins are derived (P[5] = 12 - sum, H[6] = N/2-10 - sum).
    // @pre  sums: pentagons(ni) == 12 && hexagons(ni) == N/2 - 10
    // @post result.neighbour_indices(N) == ni
    NeighbourIndices neighbour_indices(int N) const;
    static Entry with_neighbour_indices(Entry e, const NeighbourIndices& ni);
    // @pre  ascending: rspi_zero_based.size() == 12 && is_sorted(rspi_zero_based) && rspi_zero_based[0] >= 0
    static Entry with_rspi(Entry e, const vector<int>& rspi_zero_based);
    static Entry with_point_group(Entry e, const PointGroup& pg);
  };

  // What Entry::from_dual produced, or why it could not.  Only Ok leaves
  // `entry` meaningful.
  struct EntryResult {
    enum class Code {
      Ok,
      NoRegularSpiral,    // every start needs a jump: the isomer has no ring
                          // spiral of the kind the database keys on
      NotRepresentable,   // a value no record column can hold (a face position
                          // or NMR orbit above 255, a Hamilton count above I7,
                          // more than 3 NMR orbit sizes)
      UnknownPointGroup,  // Symmetry did not classify it as one of the 28
                          // fullerene point groups
    };
    Code   code = Code::Ok;
    Entry  entry{};
    string why;                      // the reason; nonempty iff code != Ok
    bool   ok() const { return code == Code::Ok; }
    static const char* code_name(Code c);
    const char* name() const { return code_name(code); }
    static EntryResult error(Code c, string why) { return {c, Entry{}, std::move(why)}; }
  };

  // The name of one field, and the fields of a set, for diagnostics.
  static const char* field_name(Field f);
  static string      field_names(Field fields);      // "group NMR"
  // The single fields of a set, in bit order -- the iteration word.
  static vector<Field> fields_of(Field fields);

  vector<Entry> entries;
  map<RSPI, int> RSPIindex;      // RSPI -> index into entries

  // The entry with this RSPI, or nullptr; and its 0-based position, or -1.
  const Entry* find(const RSPI& rspi) const {
    const int i = index_of(rspi);
    return i < 0? nullptr : &entries[i];
  }
  int index_of(const RSPI& rspi) const {
    auto it = RSPIindex.find(rspi);
    return it == RSPIindex.end()? -1 : it->second;
  }
  // The columns this file carries: every one but the Hamilton count when the
  // header says it has none.
  Field stored_fields() const { return with_ncycham? Field::All : Field::All & ~Field::Ncycham; }

  static string database_path;

  // The write routines report the I/O outcome by value (a modeled outcome:
  // the caller chose the destination); a data error -- a value that does not
  // fit its field, an entry the layout cannot hold -- throws
  // std::invalid_argument.  The read routines THROW on a missing or corrupt
  // database file (an environment error -- there is no valid IsomerDB to
  // return, and continuing with a sentinel caused silent downstream
  // corruption).  The database root comes from the FULLERENE_DATABASE_PATH
  // CMake variable (via config.hh) and can be overridden at runtime through
  // IsomerDB::database_path.
  enum class WriteStatus { Ok, CannotOpen, WriteError };
  WriteStatus writeBinary(const string& filename) const;   // @throws invalid_argument: N does not fit the 8-bit header
  WriteStatus writeCSV(const string& filename) const;

  // @throws std::runtime_error when `filename` cannot be opened, is
  //         truncated, or its entry count disagrees with its size.
  static IsomerDB readBinary(const string& filename);
  // Random access to one isomer (1-based index, matching the database
  // convention) without loading the whole file.
  // @throws std::runtime_error on a missing or truncated database file.
  // @throws std::out_of_range when isomer is outside [1, Nisomers].
  static Entry getIsomer(int N, int isomer, bool IPR=false);
  static FullereneGraph makeIsomer(int N, const Entry& e);

  // Read DB in binary format.
  // @throws as readBinary(filename).
  static IsomerDB readBinary(int N=20, bool IPR=false, string extension = "");
  static string   binaryfilename(int N, bool IPR=false, string extension = "");

  // Text database.  Layout (isomer.f formats 1003/1004/1007/1008/1009):
  //   header  I3,2I1                            N, IP (1: IPR-only file), IH (1: with Hamilton-cycle counts)
  //   record  A3,12I3,5I2,6I2,I2,I1,F7.5[,I7],6I3   IP=0: group, RSPI, PNI k=0..4, HNI k=0..5,
  //                                                  NeHOMO, NedgeHOMO, HLgap[, ncycham], INMR
  //           A3,12I3,3I2,I2,I1,F7.5[,I7],6I3       IP=1: as above with only HNI k=3..5 stored
  //                                                  (an IPR isomer has PNI=(12,0,0,0,0), HNI k<3 = 0)
  // readPDB(filename) also accepts the VERBOSE isomer listing printed by the
  // Fortran program (spiral.f/isomer.f formats 607/608, header I5,2I2, either
  // trimmed to header + records as the legacy All/c110-c120 files were, or
  // the raw program output with its banner up to the ' Isomer List Start'
  // marker): the extra derived columns (PNI(5), Np, HNI(6), sigma_h,
  // closed/open, isomer number) are checked against the stored ones and
  // dropped, so the result is the same IsomerDB the compact file would give.
  // Preamble and trailer lines (table header, summary statistics) are
  // ignored.
  //
  // Every record is held to Entry::invariant_violation (@ref
  // isomer-record-invariants), the same predicate writePDB refuses to write
  // a violation of.
  //
  // @throws std::runtime_error on a missing file, a malformed header, any
  //         record that does not parse column-exactly or violates an
  //         invariant, a repeated RSPI, or (the (N, IPR) overload, standard
  //         files only) a header disagreeing with the request or an entry
  //         count differing from number_isomers().
  static IsomerDB readPDB(const string& filename);
  static IsomerDB readPDB(int N=20, bool IPR=false, string extension = "");
  static string   PDBfilename(int N, Corpus corpus, string extension = "");
  static string   PDBfilename(int N, bool IPR=false, string extension = "") { return PDBfilename(N, IPR? Corpus::IPR : Corpus::All, extension); }
  // The sizes whose text file exists under database_path for this corpus.
  static vector<int> available_sizes(Corpus corpus);

  // Write the compact text format -- byte-identical to what the Fortran
  // writer (util.f CompressDatabase) produces, so readPDB/writePDB round-trip
  // every compact database file unchanged and convert the verbose ones.
  // Written to filename+".part" and renamed on success: no partial files.
  // @throws std::invalid_argument when a value does not fit its field width
  //         (Fortran would write asterisks; we refuse to lose data), or when
  //         IPR is set but an entry is not an IPR isomer (the IPR layout
  //         would drop its PNI/HNI columns).
  WriteStatus writePDB(const string& filename) const;

  static int64_t        number_isomers(int N, const string& sym="Any", bool IPR=false);
  static vector<string> symmetries(int N, bool IPR=false);

  static vector<size_t> Nisomers_data[2];
  static vector< vector<string> > symmetries_data[2];
  static vector< vector<size_t> > symmetry_count_data[2];

  // Field is a set: the usual algebra, so a caller can say
  // stored_fields() & ~Field::RSPI without casting.
  friend constexpr Field operator|(Field a, Field b) { return Field(uint32_t(a) | uint32_t(b)); }
  friend constexpr Field operator&(Field a, Field b) { return Field(uint32_t(a) & uint32_t(b)); }
  friend constexpr Field operator~(Field a)          { return Field(~uint32_t(a) & uint32_t(Field::All)); }
  friend constexpr bool  operator!(Field a)          { return uint32_t(a) == 0; }

  IsomerDB(int N=-1, bool IPR = false, bool IH=false,
	   vector<Entry> entries=vector<Entry>()) :
    N(N), Nisomers((int)entries.size()), IPR(IPR), with_ncycham(IH),
    entries(entries) { }

};

#endif
