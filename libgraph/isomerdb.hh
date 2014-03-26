#ifndef ISOMERDB_HH
# define ISOMERDB_HH
#include <vector>
#include <map>
#include <string>
#include <assert.h>
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

  bool writeBinary(const string filename)
  {
    FILE *f = fopen(filename.c_str(),"wb");
    if(!f){
      cerr << "Couldn't open database file " << filename << " for writing: " << strerror(errno) << ".\n";
      return false;
    }
    u_int16_t header = N | (IPR<<8) | (with_ncycham<<9);
    fwrite(&header,2,1,f);

    fwrite(&Nisomers,4,1,f);
    fwrite(&entries[0],sizeof(Entry),Nisomers,f);

    return true;
  }

  static string CSVarray(u_int8_t *v, int l){ 
    string s(to_string(vector<int>(v,v+l)));
    return s.substr(1,s.size()-2);
  }

  bool writeCSV(const string filename)
  {
    ofstream f(filename.c_str());
    if(!f){
      cerr << "Couldn't open database file " << filename << " for writing: " << strerror(errno) << ".\n";
      return false;
    }
    
    f << "\"SYMGROUP\", \"RSPI\",\"PNI\",\"HNI\",\"NEHOMO\",\"NEDGEHOMO\",\"HLGAP\",\"NCYCHAM\",\"INMR\"\r\n";
    for(int i=0;i<Nisomers;i++){
      Entry e(entries[i]);
      f << "\""<<string(e.group,3) << "\",\"" << CSVarray(e.RSPI,12) << "\",\"" << CSVarray(e.PNI,5) << "\",\"" << CSVarray(e.HNI,6) << "\","
	<< int(e.NeHOMO) << "," << int(e.NedgeHOMO) << "," << e.HLgap << "," << e.ncycham << ",\"" << CSVarray(e.INMR,6) << "\"\r\n";
    }
    return true;
  }

  static IsomerDB readBinary(const string filename){
    FILE *f = fopen(filename.c_str(),"rb");
    if(!f){
      cerr << "Couldn't open database file " << filename << " for reading: " << strerror(errno) << ".\n";
      return IsomerDB(-1);
    }
    u_int16_t header;
    size_t nread1 = fread(&header,2,1,f);
    
    IsomerDB DB(header & 0xff, header >> 8 & 1, header >> 9 & 1);
    size_t nread2 = fread(&DB.Nisomers,4,1,f);
    DB.entries.resize(DB.Nisomers);
    size_t nread3 = fread(&DB.entries[0],sizeof(Entry),DB.Nisomers,f);

    if(!nread1 || !nread2 || !nread3) return IsomerDB(-1);
    return DB;
  }
 
  static Entry getIsomer(int N, int isomer, bool IPR=false){
    string filename = FULLERENE_DATABASE"/binary/c"+pad_string(to_string(N),3)+(IPR?"IPR":"all")+".bin";
    FILE *f = fopen(filename.c_str(),"rb");
    if(!f){
      cerr << "Couldn't open database file " << filename << " for reading: " << strerror(errno) << ".\n";
      abort();
    }
    u_int16_t header;
    u_int32_t Nisomers;
    Entry e;

    if(fread(&header,2,1,f) != 1){ cerr << "Read error from database file " << filename << ": " << strerror(errno) << ".\n"; }
    if(fread(&Nisomers,4,1,f) != 1){ cerr << "Read error from database file " << filename << ": " << strerror(errno) << ".\n"; }
    assert(isomer <= Nisomers);
    fseek(f,(isomer-1)*sizeof(Entry),SEEK_CUR);
    
    if(fread(&e,sizeof(Entry),1,f) != 1){ cerr << "Read error from database file " << filename << ": " << strerror(errno) << ".\n"; }

    return e;
  }

  static FullereneGraph makeIsomer(int N, const Entry& e)
  {
    vector<int> RSPI(e.RSPI,e.RSPI+12);
    for(int i=0;i<12;i++) RSPI[i]--;
    //    cout << "creating C"<<N<< " from spiral indices " << RSPI << endl;
    return FullereneGraph(N,RSPI);
  }

  static IsomerDB readBinary(int N=20, bool IPR=false) {
    string filename = FULLERENE_DATABASE"/binary/c"+pad_string(to_string(N),3,'0')+(IPR?"IPR":"all")+string(".bin");
    return readBinary(filename);
  }

  static IsomerDB readPDB(int N=20, bool IPR=false) {
    string filename;
    if(IPR)
      filename= "database/IPR/c"+pad_string(to_string(N),3,'0')+"IPR.database";      
    else 
      filename= "database/All/c"+pad_string(to_string(N),3,'0')+"all.database";

    ifstream dbfile(filename.c_str());
    if(!dbfile){
      cerr << "Couldn't open database file " << filename << " for reading. (error: " << strerror(errno) << ")\n";
      return IsomerDB(-1);
    }
    
    string line;
    int Nread,IP,IH, pos=0;
    getline(dbfile,line);
    
    Nread = fortran_readI(line,pos,3);
    IP = fortran_readI(line,pos,1);
    IH = fortran_readI(line,pos,1);

    IsomerDB DB(Nread,IP,IH);
    
  /* 
     IP=0,IH=1:  1004 Format(A3,12I3,5I2,6I2,I2,I1,F7.5,I7,6I3)
     IP=0,IH=0:  1007 Format(A3,12I3,5I2,6I2,I2,I1,F7.5,6I3)
     IP=1,IH=1:  1008 Format(A3,12I3,3I2,I2,I1,F7.5,I7,6I3)
     IP=1,IH=0:  1009 Format(A3,12I3,3I2,I2,I1,F7.5,6I3)
  */
    while(getline(dbfile,line)){
      int pos = 0;
      Entry e;
      fortran_readA(e.group,line,pos,3);
      for(int i=0;i<12;i++) e.RSPI[i] = fortran_readI(line,pos,3);
      if(IP==0){
	for(int i=0;i<5;i++) e.PNI[i] = fortran_readI(line,pos,2);
	for(int i=0;i<3;i++) e.HNI[i] = 0;
	for(int i=3;i<6;i++) e.HNI[i] = fortran_readI(line,pos,2);
	// for(int i=0;i<6;i++) e.HNI[i] = fortran_readI(line,pos,2); // Discrepancy between util.f and db files
      } else {
	for(int i=0;i<3;i++) e.HNI[i] = 0;
	for(int i=3;i<6;i++) e.HNI[i] = fortran_readI(line,pos,2);
      }
      e.NeHOMO = fortran_readI(line,pos,2);
      //      e.NedgeHOMO = fortran_readI(line,pos,1); // Discrepancy between util.f and db files
      e.NedgeHOMO = fortran_readI(line,pos,2);
      fortran_readI(line,pos,2); // Mysterious column in db files
      e.HLgap = fortran_readF(line,pos,10);
      if(IH==1) e.ncycham = fortran_readI(line,pos,7);
      for(int i=0;i<6;i++) e.INMR[i] = fortran_readI(line,pos,3);
      DB.entries.push_back(e);
    }
    
    DB.Nisomers = DB.entries.size();
    for(int i=0;i<DB.Nisomers;i++){
      const Entry &e(DB.entries[i]);
      DB.RSPIindex[vector<int>(e.RSPI,e.RSPI+12)] = i;
    }

    return DB;
  }

  static size_t number_isomers(int N, const string& sym="Any"){ 
    int Nindex = (N-20)/2;
    if(Nindex >= Nisomers_data.size()) return 0;

    if(sym == "Any") return Nisomers_data[Nindex]; 

    if(sym == "Nontrivial"){
      size_t sum = 0;
      for(int i=0;i<symmetries_data[Nindex].size();i++) 
	if(symmetries_data[Nindex][i] != " C1") 
	  sum += symmetry_count_data[Nindex][i];
      return sum;
    } else {
      for(int i=0;i<symmetries_data[Nindex].size();i++) 
	if(symmetries_data[Nindex][i] == sym) 
	  return symmetry_count_data[Nindex][i];
    }
    return 0;
  }
  static vector<string> symmetries(int N){ return symmetries_data[(N-20)/2]; }

  static vector<size_t> Nisomers_data;
  static vector< vector<string> > symmetries_data;
  static vector< vector<size_t> > symmetry_count_data;

  IsomerDB(int N=-1, bool IPR = false, bool IH=false, 
	   vector<Entry> entries=vector<Entry>()) : 
    N(N), IPR(IPR), with_ncycham(IH), entries(entries) { }

};

#endif
