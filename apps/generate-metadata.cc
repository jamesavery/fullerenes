#include <contrib/buckygen-wrapper.hh>
#include <libgraph/triangulation.hh>
#include <libgraph/symmetry.hh>

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
  int Nmax = strtol(av[1],0,0), IPR = strtol(av[2],0,0), only_nontrivial = strtol(av[3],0,0);
  int first = IPR?60:20;
  vector<size_t> Nisomers((Nmax-first)/2+1);
  vector< map<string,size_t> > symmetry_count((Nmax-first)/2+1);
  vector< vector<string> >      symmetries((Nmax-first)/2+1);
  vector< vector<size_t> >   symmetry_count_data((Nmax-first)/2+1);

  for(int N=first;N<=Nmax;N+=2){
    int Nindex = (N-first)/2;
    if(N==22){
      Nisomers[Nindex] = 0;
      continue;
    }

    BuckyGen::buckygen_queue BQ = BuckyGen::start(N,IPR,only_nontrivial);
    Triangulation g;

    size_t cnt = 0;
    while(BuckyGen::next_fullerene(BQ,g)){
      Triangulation::jumplist_t jumps;
      vector<int> spiral(g.N/2+2);
      cnt++;
      //      g.update_from_neighbours();
      //      cout << "Calling spiral on graph " << g << "\n";
      if(g.get_spiral(spiral,jumps,false,false,false)){
	Symmetry S(spiral);
	symmetry_count[Nindex][S.point_group().to_string()]++;
	//	cout << "spiral = " << spiral << ";\n";
	//	   << "jumps  = " << jumps  << ";\n";
      } else cout << "spiral failed for " << g << endl;
    }
    BuckyGen::stop(BQ);

    Nisomers[Nindex] = cnt;
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
