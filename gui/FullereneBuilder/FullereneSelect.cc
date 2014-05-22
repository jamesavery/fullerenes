#include "FullereneSelect.hh"

FullereneSelect::FullereneSelect(int N, bool IPR, string extension) {
  // Is there a binary DB?
  has_db = true;
  printf("FullereneSelect(%d,%d,%s)\n",N,IPR,extension.c_str());
  db = IsomerDB::readBinary(N,IPR,extension);
  if(db.N<0) db = IsomerDB::readPDB(N,IPR,extension);
  if(db.N<0) { db.N = N; db.IPR = IPR; has_db = false; }
}



vector<Entry> FullereneSelect::get_fullerenes(int max_results, int iso_start, int iso_end, string sym_filter, int oc_filter, double gap_min, double gap_max) const
{
  iso_end = min(iso_end,iso_start+max_results-1);
  vector<Entry> results;
  if(has_db){
    for(int i=iso_start; i<=iso_end && i <= db.entries.size(); i++){
      const IsomerDB::Entry &e = db.entries[i];
      if(sym_filter != "Any" && string(e.group,3) != sym_filter){
	//	printf("sym-filter '%s' != '%s'\n",string(e.group,3).c_str(),sym_filter.c_str());
	continue;
      }
      //      if((oc_filter & 2) && // How to check open/closed?
      if(e.HLgap < gap_min || e.HLgap > gap_max){
	//	printf("gap-filter %g %g %g\n",gap_min,e.HLgap,gap_max);
	continue;
      }
      
      results.push_back(Entry(i,e));
    }

  } else {
    // Generate the lot - store graphs, but don't do RSPI, perhaps?
  }
 return results;

}

vector<double> FullereneSelect::get_coordinates(const FullereneGraph &g, int opt_method, double tolerance) const
{
  vector<coord3d> geom0 = g.zero_order_geometry();
  vector<coord3d> geom  = g.optimized_geometry(geom0,opt_method,tolerance);
  vector<double>  coords(db.N*3);
  for(int i=0;i<db.N;i++){
    coords[i*3+0] = geom[i].x[0];
    coords[i*3+1] = geom[i].x[1];
    coords[i*3+2] = geom[i].x[2];
  }
  return coords;
}  

vector<double> FullereneSelect::get_coordinates(int isomer, int opt_method, double tolerance) const
{
  FullereneGraph g;
  if(has_db){
    g = IsomerDB::makeIsomer(db.N,db.entries[isomer]);
    g.layout2d = g.tutte_layout();
  } else {
    // get generated graph
  }
  return get_coordinates(g,opt_method,tolerance);
}

Triangulation FullereneSelect::get_triangulation(const Entry& e) const
{
  int N = db.N;
  vector<int> spiral(N/2+2,6);
  for(int i=0;i<12;i++) spiral[e.RSPI[i]-1] = 5;
  return Triangulation(spiral);
}


