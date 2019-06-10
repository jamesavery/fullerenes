#ifndef FULLERENEGUI_HH
# define FULLERENEGUI_HH

#include <libgraph/isomerdb.hh>
#include <libgraph/triangulation.hh>
#include <contrib/buckygen-wrapper.hh>
#include <vector>


using namespace std;

struct Entry {
  int isomer;
  string group;		
  vector<int> RSPI;		
  vector<int> PNI;		
  vector<int> HNI;		
  int NeHOMO;		
  int NedgeHOMO;		
  float HLgap;		
  int ncycham;		
  vector<int> INMR;

  Entry() {}
  Entry(int isomer,const IsomerDB::Entry& e) : isomer(isomer),
    group(e.group,3), RSPI(e.RSPI,e.RSPI+12), PNI(e.PNI,e.PNI+5), HNI(e.HNI,e.HNI+6),
    NeHOMO(e.NeHOMO), NedgeHOMO(e.NedgeHOMO), HLgap(e.HLgap), ncycham(e.ncycham),
    INMR(e.INMR,e.INMR+6) {}
};

class FullereneSelect {
public:
  IsomerDB db;
  bool has_db;

  FullereneSelect(int N, bool IPR, string extension="");

  vector<Entry> get_fullerenes(int max_results, int iso_from, int iso_to, string sym_filter="Any", int oc_filter=0, double gap_min=-1e6, double gap_max=1e6) const;

  Triangulation get_triangulation(const Entry& e) const;
  vector<double> get_coordinates(const FullereneGraph &g, int opt_method, double tolerance) const;
  vector<double> get_coordinates(int isomer, int opt_method, double tolerance) const;


  
};

#endif
