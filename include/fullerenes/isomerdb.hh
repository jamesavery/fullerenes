#ifndef ISOMERDB_HH
# define ISOMERDB_HH
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <errno.h>

#include "fullerenegraph.hh"

using namespace std;

class IsomerDB { 
public:
  int N, Nisomers;
  bool IPR, with_ncycham;

  struct Entry {
    char group[3];		
    u_int8_t RSPI[12];		
    u_int8_t PNI[5];		
    u_int8_t HNI[6];		
    u_int8_t NeHOMO;		
    u_int8_t NedgeHOMO;		
    float HLgap;		
    int ncycham;		
    u_int8_t INMR[6];		
  };

  vector<Entry> entries;
  map<vector<int>, int> RSPIindex;

  static string database_path;

  static void fortran_readA(char *result, string s, int& pos, int len){
    for(int i=0;i<len;i++) result[i] = s[pos+i];
    pos += len;
  }

  static int fortran_readI(string s, int& pos, int len){
    int result = strtol(s.substr(pos,len).c_str(),0,0);
    pos += len;
    return result;
  }

  static double fortran_readF(string s, int& pos, int len){
    double result = strtod(s.substr(pos,len).c_str(),0);
    pos += len;
    return result;
  }

  static string CSVarray(u_int8_t *v, int l){ 
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
  // Read DB in Peter's ASCII text format.
  // @throws std::runtime_error on a missing database file or a malformed
  //         header line.
  static IsomerDB readPDB(int N=20, bool IPR=false, string extension = "");

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
