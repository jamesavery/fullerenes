#include "libgraph/fullerenegraph.hh"
//#include "libgraph/polyhedron.hh"
#include "libgraph/triangulation.hh"


// triangulations: we get the permutation of vertices. done.
// cubic graphs: get dual (how are the new vertex numbers derived from the old?), order faces, translate back (how are indices inherited?)
// leap-frog: (permutation of leap frog?), identify "


// yes, I know this is inefficient, but the vectors are short
int intersection_3(vector<int> &a, vector<int> &b, vector<int> &c){
  for (int i=0; i<a.size(); i++){
    for (int j=0; j<b.size(); j++){
      if (a[i] != b[j]) continue;
      for (int k=0; k<c.size(); k++){
        if (a[i] == c[k]) return a[i];
      }
    }
  }
  assert(true);
  return -1; // not good
}


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

  // get triangulation with face numbers (wrt cub)
  Triangulation T2 (FG.dual_graph());

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
  cout << "triangles, relative to trig internal vertex numbers: " << tri_faces_orig << endl;
  vector<tri_t> tri_faces_permuted(tri_faces_orig);
  for(int i=0; i<tri_faces_orig.size(); i++){
    for(int j=0; j<3; j++){
      tri_faces_permuted[i][j] = permutation2_inv[tri_faces_permuted[i][j]];
    }
  }
  // sort triples, first the vertices per face, than the faces
  for(tri_t& t: tri_faces_permuted){
    if(t[0] > t[1]) swap(t[0], t[1]);
    if(t[1] > t[2]) swap(t[1], t[2]);
    if(t[0] > t[1]) swap(t[0], t[1]);
  }
  sort(tri_faces_permuted.begin(), tri_faces_permuted.end()); 
  cout <<  "triangles, relative to canon vertex numbers: " << tri_faces_permuted << endl;
  
  // permute back
  for(int i=0; i<tri_faces_orig.size(); i++){
    for(int j=0; j<3; j++){
      tri_faces_permuted[i][j] = permutation2[tri_faces_permuted[i][j]];
    }
  }
  cout << "triangles, in terms of internal numbers, sorted according to canon numbers: " << tri_faces_permuted << endl;

  // and then find the vertex in cub, which is part of faces {0,1,2},{0,1,5}, ...
  vector<face_t> faces = FG.compute_faces_oriented();
  cout << "faces " << faces << endl;
  vector<int> vertex_numbers_inv(FG.N);
  for(int v=0; v<FG.N; v++){
    face_t f0 = faces[tri_faces_permuted[v][0]];
    face_t f1 = faces[tri_faces_permuted[v][1]];
    face_t f2 = faces[tri_faces_permuted[v][2]];
    vertex_numbers_inv[v] = intersection_3(f0, f1, f2);
  }
  cout << "permutation of cub vertex numbers inv: " << vertex_numbers_inv << endl;
  vector<int> vertex_numbers(FG.N);
  for(int i=0; i< vertex_numbers.size(); i++) vertex_numbers[vertex_numbers_inv[i]] = i;
  cout << "permutation of cub vertex numbers (ie, replace v by vertex_numbers[v], to get numbered vertices): " << vertex_numbers << endl;


  return 0;
}
