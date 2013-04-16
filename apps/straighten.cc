#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include "libgraph/cubicgraph.hh"
#include "libgraph/fullerenegraph.hh"
#include "libgraph/polyhedron.hh"

struct approx_sort : public std::binary_function<coord3d, coord3d, bool>
{
  bool operator()(const coord3d &x, const coord3d &y) const
  {   
    coord3d d(x-y);
    for(int i=0;i<3;i++) d.x[i] = round(d.x[i]*pow(10,2))/pow(10,2);

    return d.x[0] < 0 || ((d.x[0] == 0 && d.x[1] < 0) || (d.x[1] == 0 && d.x[2] < 0));
  }
};

vector<coord3d> readxyz(const string filename)
{
  ifstream xyzfile(filename.c_str());
  string line,comment,dummy;
  int N, linenumber = 0,k=0;


  if(xyzfile.fail()){
    cerr << "Couldn't open file " << filename << " for reading.\n";
    abort();
  }
  cout << "Reading " << filename << endl;
  set<coord3d,approx_sort> coordinate_set;

  while(getline(xyzfile,line)){  
    stringstream l(line);
    ++linenumber;
    //    cout << linenumber << ": " << line << endl;
    if(linenumber == 1) {
      l >> N;
    }
    if(linenumber == 2) comment = line;
    if(linenumber > 2) {
      coord3d x;
      l >> dummy;
      for(int i=0;i<3;i++)
	l >> x.x[i];
      if(l.fail()){ cerr << "Malformed line: " << line << endl;  continue; }      
      coordinate_set.insert(x);
      k++;
    }    
  }
  xyzfile.close();

  vector<coord3d> coordinates(coordinate_set.begin(),coordinate_set.end());    
  cout << "coordinates = " << coordinates << endl;
  cout << "coordinates.size() = " << coordinates.size() << endl;
  assert(coordinates.size() == N);


  return coordinates;
}

double diameter(const vector<coord3d> &xs)
{
  double dmax = -INFINITY;
  for(int i=0;i<xs.size();i++)
    for(int j=i+1;j<xs.size();j++){
      double d = (xs[i]-xs[j]).norm();
      if(d>dmax) dmax = d;
    }
  return dmax;
}

coord3d centroid(const vector<coord3d> &xs)
{
  coord3d centroid;
  for(int i=0;i<xs.size();i++)
    centroid += xs[i];
  centroid *= 1.0/double(xs.size());
  return centroid;
}

int main(int ac, char **av)
{
  if(ac<4){
    cerr << "Syntax: " << av[0] << " <src.xyz>  <tgt.xyz> <out.xyz>\n";
    return -1;
  }

  vector<coord3d> src(readxyz(av[1])), tgt(readxyz(av[2]));

  assert(src.size() == tgt.size());

  coord3d c_src(centroid(src)), c_tgt(centroid(tgt));
  double d_src(diameter(src)), d_tgt(diameter(tgt));

  cout << "Translating from centroid " << c_src << " to " << c_tgt << endl;
  cout << "and scaling from diameter " << d_src << " to " << d_tgt << endl;
    
  vector<coord3d> result(src.size());

  for(int i=0;i<src.size();i++)
    result[i] = (src[i] - c_src)*(d_tgt/d_src) + c_tgt;

  ofstream outfile(av[3]);
  outfile << " " << src.size() << endl;
  outfile << "Auto-transformed coordinates\n";
  for(int i=0;i<src.size();i++)
    outfile << " C " << result[i][0] << " " << result[i][1] << " " << result[i][2] << endl;

  return 0;
}


