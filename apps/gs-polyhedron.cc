#include "libgraph/fullerenegraph.hh"
#include "libgraph/polyhedron.hh"
#include "libgraph/isomerdb.hh"
#include "libgraph/triangulation.hh"
#include "libgraph/fullerenegraph.hh"
#include "libgraph/symmetry.hh"

// all lists are CCW

Polyhedron polyhedron_from_spiral(const general_spiral& GS, bool tri=false)
{
  Triangulation g(GS.spiral,GS.jumps);
  int max_face = 0;
  for(auto f: GS.spiral) max_face = max(f,max_face);

  PlanarGraph cubic(g.dual_graph());
  cubic.layout2d = cubic.tutte_layout(-1,-1,-1,max_face);
  vector<coord3d>  points = cubic.zero_order_geometry();
  vector<face_t> faces = cubic.compute_faces();

  Polyhedron P(cubic,points,max_face,faces);
  P.optimize();

  return P;
}


Polyhedron LFPolyhedron(const Polyhedron& P)
{
  Triangulation LF = P.leapfrog_dual();
  vector<coord3d> points = P.points;
  points.resize(LF.N);

  for(int i=0;i<P.faces.size();i++){
    points[P.N+i] = P.faces[i].centroid(P.points);
  }

  return Polyhedron(LF,points);
}

string jumps_to_string(const jumplist_t &jumps)
{
  string s="";
  for(const auto &j: jumps)
    s += to_string(j.first) +"," + to_string(j.second) + ",";
  s.pop_back();

  return s;
}

string spiral_to_string(const vector<int>& spiral)
{
  string s="";
  for(int i=0;i<spiral.size();i++)
    s += to_string(spiral[i]) + (i+1<spiral.size()? ",":"");
  return s;
}

vector<int> spiral_to_rspi(vector<int> spiral)
{
  vector<int> rspi(12);
  for(int i=0,j=0;i<spiral.size();i++) if(spiral[i] == 5) rspi[j++] = i;
  return rspi;
}


string spiral_to_rspi_string(const vector<int>& spiral)
{
  return spiral_to_string(spiral_to_rspi(spiral)+1);
}

int main(int ac, char **av)
{
  general_spiral GSin, GSout, CSout;

  int i=1;
  for(;i<ac && av[i][0] != '-';i+=2){
    GSin.jumps.push_back(make_pair((int)strtol(av[i],0,0)-1,(int)strtol(av[i+1],0,0)));
  }
  if(i>=ac){
    fprintf(stderr,"Syntax: %s <jumps> '-' <spiral>\n"
	    " where <jumps> ::= n_1 k_1 ... n_m k_m and <spiral> ::= f_1 ... f_N\n\n",
	    av[0]);
    return -1;
  }
  for(i++;i<ac; i++)
    GSin.spiral.push_back(strtol(av[i],0,0));




  Triangulation triangulation(GSin.spiral,GSin.jumps);
  PlanarGraph   cubic(triangulation.dual_graph());

  string basename = "C"+to_string(cubic.N)+"-"+spiral_to_string(GSin.spiral);  
  ofstream output("output/"+basename+".m");
  output << "GSin  = " << GSin << ";\n";
  

  
  triangulation.get_spiral(GSout.spiral,GSout.jumps);
  bool unwind_success = triangulation.get_spiral(CSout.spiral,CSout.jumps,false,true);
  
  if(unwind_success){
    output << "GSout = " << GSout << ";\n";
    output << "CSout = " << CSout << ";\n";    
  } else {
    cerr << "Spiral unwinding failed.\n";
    return -1;
  }

  output << "dg = " << triangulation << ";\n"
	 << "g  = " << cubic << ";\n";
  output.flush();

  Polyhedron P(polyhedron_from_spiral(GSin));

  {
    ofstream mol2(("output/"+basename+"-P.mol2").c_str());
    mol2 << P.to_mol2();
    mol2.close();

    ofstream pov(("output/"+basename+"-P.pov").c_str());
    pov << P.to_povray();
    pov.close();
  }
  
  output.close();
  return 0;
}
