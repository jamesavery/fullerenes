#include "fullerenes/spiral.hh"
#include "fullerenes/fullerenegraph.hh"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/isomerdb.hh"
#include "fullerenes/unfold.hh"
#include "fullerenes/eisenstein.hh"

#include <teem/nrrd.h>

struct volume {
  coord3d X0, X1;
  double dx;
  int nx[3];

  vector<float> data;

  void to_raw(string filename)
  {
    FILE *f = fopen(filename.c_str(),"w");
    fwrite(&X0.x[0],sizeof(double),3,f);
    fwrite(&X1.x[0],sizeof(double),3,f);
    fwrite(&dx,sizeof(double),1,f);
    fwrite(&nx[0],sizeof(int),3,f);
    fwrite(&data[0],sizeof(float),nx[0]*nx[1]*nx[2],f);
    fclose(f);
  }
};

// Fake density -- todo: Read real density from HDF5
float density_eval(const Polyhedron &P, const coord3d &X)
{
  double sum = 0;
  for(auto Xp: P.points){
    double r = (Xp-X).norm();
    sum += exp(-1.7*r);
  }
  return sum;
}

volume density(const Polyhedron &P, float dx)
{
  coord3d X0 = {INFINITY, INFINITY, INFINITY}, X1 = {-INFINITY, -INFINITY, -INFINITY};
  for(auto p: P.points){
    for(int i=0;i<3;i++){ X0[i] = min(X0[i],2*p[i]); X1[i] = max(X1[i],2*p[i]); }
  }

  int nx = ceil((X1[0]-X0[0])/dx);
  int ny = ceil((X1[1]-X0[1])/dx);
  int nz = ceil((X1[2]-X0[2])/dx);
  
  vector<float> data(nx*ny*nz);

  float x,y,z;
  int i,j,k;
  
#pragma omp parallel for collapse(3)
  for(i=0, x=X0[0];i<nx;i++,x+=dx)
    for(j=0, y=X0[1];j<ny;j++, y+=dx)
      for(k=0, z=X0[2];k<nz;k++, z+= dx){
	//	printf("%g,%g,%g\n",x,y,z);
	data[i*ny*nz + j*nz + k] = density_eval(P,{x,y,z});
      }

  return volume{X0, X1, dx, {nx,ny,nz}, data};
}



vector< pair<Eisenstein,node_t> > outline{
  {{0,14},  1 }, {{0,16},  12}, {{2,16},  5 }, {{4,16},   13},
  {{6,14},  3 }, {{8,12},  14}, {{10,12}, 2 }, {{8,14},   14},{{8,16},3},
  {{10,16}, 12}, {{12,14}, 4 }, {{12,16}, 12}, {{14,16},  5 },
  {{16,14}, 13}, {{16,12}, 9 }, {{16,10}, 15},		   
  {{14,10}, 7 }, {{12,10}, 16}, {{12, 8}, 6 }, {{14, 8}, 16},{{16,6},7},
  {{16, 4}, 17}, {{14, 4}, 8 }, {{16, 2}, 17}, {{16, 0}, 9},
  {{14, 0}, 18}, {{12, 2}, 1 }, {{10, 4}, 19},		   
  {{10, 6}, 11}, {{10, 8}, 20}, {{8, 10}, 10}, {{8, 8}, 20},{{6, 8},11},
  {{4, 10}, 19}, {{4, 12}, 0 }, {{2, 12}, 19}
};

vector<coord2d> uvmap(const PlanarGraph& G)
{
  Triangulation dG = G.dual_graph();
  Unfolding U(dG);

  // Grid range
  int minx=INT_MAX, miny=INT_MAX, maxx=INT_MIN, maxy=INT_MIN;
  int x,y;
  for(auto o: U.outline){
    tie(x,y) = o.first;
    minx = min(minx,x); maxx = max(maxx,x);
    miny = min(miny,y); maxy = max(maxy,y);
  }
  
  // dedge -> triangle (atom)
  node_t q,r,s,t;
  Eisenstein xq, xr, xs, xt;
  for(auto kv: U.arc_coords){
    tie(q,r)   = kv.first;
    tie(xq,xr) = kv.second;

    
  }

  vector<coord2d> surface;

  
  return surface;
}

// bool to_wavefront_obj(const Polyhedron &P, FILE *file)
// {
//   fprintf(file,"# Vertices:\n");    
//   for(auto p: P.points)
//     fprintf(file,"v %f %f %f\n",p[0],p[1],p[2]);

//   fprintf(file,"# Pentagons:\n"
// 	       "g pentagons\n");
//   for(auto f: P.faces)
//     if(f.size()==5)
//       fprintf(file,"f %d %d %d %d %d\n",f[0]+1,f[1]+1,f[2]+1,f[3]+1,f[4]+1);

//   fprintf(file,"# Hexagons:\n"
// 	       "g hexagons\n");
//   for(auto f: P.faces)
//     if(f.size()==5)
//       fprintf(file,"f %d %d %d %d %d %d\n",f[0]+1,f[1]+1,f[2]+1,f[3]+1,f[4]+1,f[5]+1);  
// }


int main(int ac, char **av)
{
  int N;
  jumplist_t jumps;
  vector<int> RSPI(12);
  bool from_file = false;

  if(ac==2){			// If only one argument is given, it is a filename to read initial polyhedron from
    from_file = true;
    N = 0;
    
  } else if(ac==3) {            // Pentagon indices in quoted string
    N = strtol(av[1],0,0);
    int n = strlen(av[2]);
    for(int i=0;i<n;i++) if(av[2][i] == ',') av[2][i] = ' ';
    istringstream iss(av[2]);
    for(int i=0;i<12;i++){ iss >> RSPI[i]; RSPI[i]--; }
    cout << "RSPI = " << RSPI << endl;
  } else 
  if(ac>=14) {
    N = strtol(av[1],0,0);
    for(int i=0;i<12;i++) RSPI[i] = strtol(av[i+2],0,0)-1;
  } 


  // int isomer_number = ac>=3? strtol(av[2],0,0) : 1;
  // bool IPR = ac>=4? strtol(av[3],0,0) : false;

  // IsomerDB DB(N,IPR);

  // FullereneGraph g(IsomerDB::makeIsomer(N,DB.getIsomer(N,isomer_number,IPR)));


  vector<int> spiral(N/2+2,6);
  for(int i=0;i<12;i++) spiral[RSPI[i]] = 5;

  Polyhedron P0;
  PlanarGraph g;
  if(from_file){
    P0 = Polyhedron::from_file(av[1]);
    N = P0.N;
    g = P0;
    g.layout2d = g.tutte_layout(-1,-1,-1,8);
  } else {
    Triangulation T(spiral,jumps);
    g = T.dual_graph();
    g.layout2d = g.tutte_layout();
    P0 = Polyhedron(g,g.zero_order_geometry(),6);
  }

  string basename("polyhedron-"+to_string(N));
  Polyhedron::to_file(P0,"output/"+basename+"-P0.mol2");

  printf("P0\n");
  Polyhedron P(P0);
  printf("Optimizing P\n");  
  P.optimise();
  printf("Writing P\n");
  Polyhedron::to_file(P,"output/"+basename+".mol2");

  printf("Aligning P\n");  
  P.move_to_origin();
  P.align_with_axes();

  printf("Writing P-aligned\n");    
  Polyhedron::to_file(P,"output/"+basename+"-if.mol2");
  Polyhedron::to_file(P,"output/"+basename+"-if.obj");  

  // volume rho = density(P,0.05);
  // rho.to_raw("output/"+basename+"-density.raw");

  Polyhedron Plf = P.leapfrog_dual();
  Polyhedron::to_file(Plf,"output/"+basename+"-Plf.mol2");
  Polyhedron::to_file(Plf,"output/"+basename+"-Plf.obj");      

  vector<coord2d> texmap = uvmap(P);

  return 0;
}
