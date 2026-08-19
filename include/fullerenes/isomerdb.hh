#ifndef ISOMERDB_HH
# define ISOMERDB_HH
#include <array>
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <errno.h>

#include "fullerenegraph.hh"

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

    IsomerDB::RSPI rspi() const { IsomerDB::RSPI r; std::copy(RSPI, RSPI+12, r.begin()); return r; }
    // The spiral as the graph constructors take it (FullereneGraph(N, rspi),
    // FullereneDual(N, rspi)): 0-based positions.
    vector<int> rspi_zero_based() const { vector<int> r(RSPI, RSPI+12); for(int& x: r) x--; return r; }
    // The neighbour indices as stored: PNI k<5, HNI k<6; the dropped top
    // bins are derived (P[5] = 12 - sum, H[6] = N/2-10 - sum).
    NeighbourIndices neighbour_indices(int N) const;
    static Entry with_neighbour_indices(Entry e, const NeighbourIndices& ni);
  };

  vector<Entry> entries;
  map<RSPI, int> RSPIindex;      // RSPI -> index into entries

  // The entry with this RSPI, or nullptr.
  const Entry* find(const RSPI& rspi) const {
    auto it = RSPIindex.find(rspi);
    return it == RSPIindex.end()? nullptr : &entries[it->second];
  }

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
  // Every record is held to the invariants of a fullerene record (RSPI in
  // [1, N/2+2] ascending, PNI/HNI sums, NMR orbits summing to N, HOMO
  // occupation, gap zero exactly for an open shell, Hamilton count present
  // exactly when the header says so) -- see check_entry in isomerdb.cc.
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

  IsomerDB(int N=-1, bool IPR = false, bool IH=false,
	   vector<Entry> entries=vector<Entry>()) :
    N(N), Nisomers((int)entries.size()), IPR(IPR), with_ncycham(IH),
    entries(entries) { }

};

#endif
