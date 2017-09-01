#include "libgraph/fullerenegraph.hh"
//#include "libgraph/polyhedron.hh"
#include "libgraph/triangulation.hh"


// triangulations: we get the permutation of vertices. done.
// cubic graphs: get dual (how are the new vertex numbers derived from the old?), order faces, translate back (how are indices inherited?)
// leap-frog: (permutation of leap frog?), identify "



int testN = 28;
int testRSPI[12] = {1,2,3,4,5,7,10,12,13,14,15,16};

int main(int ac, char **av)
{
  int N;
  Triangulation::jumplist_t jumps;
  vector<int> RSPI(12);
  if(ac==2){
    N = 0;
  } else if(ac<14){
    N = testN;
    for(int i=0;i<12;i++) RSPI[i] = testRSPI[i]-1;
  }
  if(ac>=14) {
    N = strtol(av[1],0,0);
    for(int i=0;i<12;i++) RSPI[i] = strtol(av[i+2],0,0)-1;
  }
  if(ac>14){			// General RSPI: RSPI followed by jumps.
    for(int i=14;i<ac;i+=2)
      jumps.push_back(make_pair(strtol(av[i],0,0)-1,strtol(av[i+1],0,0)));
  }


  vector<int> spiral(N/2+2,6);
  for(int i=0;i<12;i++) spiral[RSPI[i]] = 5;

// vertex numbers in the triangulation
  cout << "*** trig" << endl;
  Triangulation T1(spiral,jumps);
  cout << T1 << endl;

  vector<int> S1;
  jumplist_t J1;
  vector<int> permutation1;
  T1.get_spiral(S1, J1, permutation1);
  cout << "vertex numbers of triangulation: " << permutation1 << endl;


// vertex numbers in the cubic graph
  cout << "*** cub" << endl;
  FullereneGraph FG(N, RSPI,jumps);
  FG.layout2d = FG.tutte_layout();

  unordered_map<dedge_t,int> face_numbers;
  size_t Fmax = 15;
  // get triangulation with face numbers (wrt cub)
  Triangulation T2 (FG.dual_graph(6, true, &face_numbers));
     cerr << "face_reps' = " << get_keys(face_numbers) << ";\n";
     cerr << "face_nums  = " << get_values(face_numbers) << ";\n";
        //cerr << "faces      = " << compute_faces_oriented(Fmax) << ";\n";
        //    cerr << "Processing face " << i_f << ": " << e_f << " -> " << get_face_oriented(e_f,Fmax) << ";\n";
  vector<int> S2;
  jumplist_t J2;
  vector<int> permutation2;
  // get spiral permutation (vertices in trig)
  T2.get_spiral(S2, J2, permutation2);
  cout << "permutation of vertex numbers of triangulation: " << permutation2 << endl;
  vector<int> permutation2_inv(permutation2.size());
  for(int i=0; i< permutation2.size(); i++) permutation2_inv[permutation2[i]] = i;
  cout << "inverse permutation of vertex numbers of triangulation: " << permutation2_inv << endl;

  vector<tri_t> tri_faces_orig = T2.compute_faces();
  cout << tri_faces_orig << endl;
  vector<tri_t> tri_faces_permuted(tri_faces_orig);
  for(int i=0; i<tri_faces_orig.size(); i++){
    for(int j=0; j<3; j++){
      tri_faces_permuted[i][j] = permutation2_inv[tri_faces_permuted[i][j]];
    }
  }
  cout << tri_faces_permuted << endl;
  //sort next
  
  


  // string basename("polyhedron-"+to_string(N));
  // {
  //   ofstream mol2(("output/"+basename+"-P0.mol2").c_str());
  //   mol2 << P0.to_mol2();
  //   mol2.close();
  // }


  // vector<face_t> faces(g.compute_faces(8,true));
  // output << "g = " << g << ";\n";
  // output << "coordinates0 = " << P0.points << ";\n";
  // output << "coordinates = "  << P.points << ";\n";
  // output << "faces = " << faces << ";\n"
  // 	  << "RSPI = " << RSPI << ";\n";

  // output << "P0 = " << P0 << ";\n";
  // output << "P = " << P << ";\n";

  // Polyhedron D(P.dual(6,true));
  // D.layout2d = D.tutte_layout();
  // D.faces    = D.compute_faces(3,true);
  // D.face_max = 3;
  // //   D.optimize();
  // output << "PD = " << D << ";\n";
  //
  // output.close();

  return 0;
}
