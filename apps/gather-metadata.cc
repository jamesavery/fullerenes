#include <libgraph/isomerdb.hh>

ostream& operator<<(ostream& S, const vector<string>& set)
{
  S << "{";
  for(auto e(set.begin()); e!=set.end(); ){
    S << "\"" << *e << "\"";
    e++; if(e!=set.end()) S << ",";
  }
  S << "}";
  return S;
}

int main(int ac, char **av)
{
  if(ac<3) return -1;
  int Nmax = strtol(av[1],0,0), IPR = strtol(av[2],0,0);
  
  int first = IPR?60:20;
  vector<size_t> Nisomers((Nmax-first)/2+1);
  vector< map<string,size_t> > symmetry_count((Nmax-first)/2+1);
  vector< vector<string> >      symmetries((Nmax-first)/2+1);
  vector< vector<size_t> >   symmetry_count_data((Nmax-first)/2+1);

  for(int N=first;N<=Nmax;N+=2){
    int Nindex = (N-first)/2;
    if(N==22 || (IPR && N>60 && N<70)){
      Nisomers[Nindex] = 0;
      continue;
    }

    IsomerDB db = IsomerDB::readBinary(N,IPR);
    if(db.N < 0){
      fprintf(stderr,"Can't open database for C%d\n",N);
      return -1;
    }
    Nisomers[Nindex] = db.Nisomers;
    
    for(int iso=0;iso<db.Nisomers;iso++){
      string sym_string(db.entries[iso].group,3);
      symmetry_count[Nindex][sym_string]++;
    }
    symmetries[Nindex] = get_keys(symmetry_count[Nindex]);
    symmetry_count_data[Nindex] = get_values(symmetry_count[Nindex]);
  }
  cout
    << "#include \"isomerdb.hh\"\n"
    << "vector<size_t> IsomerDB::Nisomers_data                 = " << Nisomers << ";\n\n"
    << "vector< vector<string> > IsomerDB::symmetries_data     = " << symmetries << ";\n\n"
    << "vector< vector<size_t> > IsomerDB::symmetry_count_data = " << symmetry_count_data << ";\n\n";

  return 0;
}
