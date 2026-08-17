#ifndef ISOMERDB_HH
# define ISOMERDB_HH
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <errno.h>

#include "fullerenegraph.hh"

using namespace std;

// The fullerene isomer database ("PDB" = Peter's Database: the text format
// of the legacy Fortran program by P. Schwerdtfeger et al.), one file per
// (N, IPR) holding every isomer with its point group, canonical ring-spiral
// pentagon indices, pentagon/hexagon neighbour indices, Hueckel HOMO data,
// Hamilton-cycle count and NMR pattern.  This class is the C++ reader/writer
// of that format and replaces the Fortran Printdatabase / CompressDatabase
// routines (isomer.f, util.f).
class IsomerDB {
public:
  int N, Nisomers;
  bool IPR, with_ncycham;

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
  };

  vector<Entry> entries;
  map<vector<int>, int> RSPIindex;

  static string database_path;

  // Fortran formatted-field primitives.  The text database is written and
  // read by Fortran FORMAT statements; these emulate one edit descriptor
  // (Aw / Iw / Fw.d) applied at column pos of record s, advancing pos by len.
  // Input follows Fortran BLANK='NULL': blanks anywhere in a numeric field
  // are ignored and an all-blank field reads as 0; anything else is
  // malformed.  Output right-justifies in len columns.
  // @throws std::invalid_argument on a malformed or missing input field, or
  //         on output that does not fit (where Fortran would print asterisks
  //         and silently lose the value).
  static void   fortran_readA(char *result, const string& s, int& pos, int len);
  static long   fortran_readI(const string& s, int& pos, int len);
  static double fortran_readF(const string& s, int& pos, int len);
  static string fortran_writeI(long v, int len);
  static string fortran_writeF(double v, int len, int decimals);

  static string CSVarray(const u_int8_t *v, int l){
    string s(to_string(vector<int>(v,v+l)));
    return s.substr(1,s.size()-2);
  }

  // The write routines report failure through their return value (a
  // modeled outcome: the caller chose the destination); the read routines
  // THROW on a missing or corrupt database file (an environment error --
  // there is no valid IsomerDB to return, and continuing with a sentinel
  // caused silent downstream corruption).  The database root comes from
  // the FULLERENE_DATABASE_PATH CMake variable (via config.hh) and can be
  // overridden at runtime through IsomerDB::database_path.
  bool writeBinary(const string filename) const;
  bool writeCSV(const string filename) const;

  // @throws std::runtime_error when `filename` cannot be opened, is
  //         truncated, or its entry count disagrees with its size.
  static IsomerDB readBinary(const string filename);
  // Random access to one isomer (1-based index, matching the database
  // convention) without loading the whole file.
  // @throws std::runtime_error on a missing or truncated database file.
  // @throws std::out_of_range when isomer is outside [1, Nisomers].
  static Entry getIsomer(int N, int isomer, bool IPR=false);
  static FullereneGraph makeIsomer(int N, const Entry& e);

  // Read DB in binary format.
  // @throws as readBinary(filename).
  static IsomerDB readBinary(int N=20, bool IPR=false, string extension = "");

  // Text database.  Layout (isomer.f formats 1003/1004/1007/1008/1009):
  //   header  I3,2I1                            N, IP (1: IPR-only file), IH (1: with Hamilton-cycle counts)
  //   record  A3,12I3,5I2,6I2,I2,I1,F7.5[,I7],6I3   IP=0: group, RSPI, PNI k=0..4, HNI k=0..5,
  //                                                  NeHOMO, NedgeHOMO, HLgap[, ncycham], INMR
  //           A3,12I3,3I2,I2,I1,F7.5[,I7],6I3       IP=1: as above with only HNI k=3..5 stored
  //                                                  (an IPR isomer has PNI=(12,0,0,0,0), HNI k<3 = 0)
  // readPDB(filename) also accepts the VERBOSE isomer listing printed by the
  // Fortran program (spiral.f/isomer.f formats 607/608, header I5,2I2 --
  // the format of the unconverted All/c110-c120 files): the extra derived
  // columns (PNI(5), Np, HNI(6), sigma_h, closed/open, isomer number) are
  // checked against the stored ones and dropped, so the result is the same
  // IsomerDB the compact file would give.  Trailer lines after the records
  // (the program's summary statistics) are ignored.
  //
  // @throws std::runtime_error on a missing file, a malformed header, or any
  //         record that does not parse column-exactly (a truncated line, a
  //         value that does not fit its field, an RSPI outside [1, N/2+2] or
  //         not ascending, a verbose record whose derived columns disagree).
  //         The (N, IPR) overload also throws when the file header disagrees
  //         with the request.
  static IsomerDB readPDB(const string& filename);
  static IsomerDB readPDB(int N=20, bool IPR=false, string extension = "");
  static string   PDBfilename(int N, bool IPR=false, string extension = "");

  // Write the compact text format -- byte-identical to what the Fortran
  // writer (util.f CompressDatabase) produces, so readPDB/writePDB round-trip
  // every compact database file unchanged and convert the verbose ones.
  // @throws std::invalid_argument when a value does not fit its field width
  //         (Fortran would write asterisks; we refuse to lose data).
  bool writePDB(const string& filename) const;

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
