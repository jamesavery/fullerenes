#include "cubicgraph.hh"
#include "triangulation.hh"

#include <cstdlib>
#include <iostream>
#include <vector>
#include <fstream>


CubicGraph::CubicGraph(FILE *file){
  char line[0x300];
  while(!feof(file)){
    node_t n, ns[3];
    double x, y;
    char *p = fgets(line,0x2ff,file);
    if(!p){
      if(feof(file)) continue;
      else {
        fprintf(stderr,"File read error.\n");
        abort();
      }
    }

    int count = sscanf(line,"%d %lf %lf %d %d %d",&n,&x,&y,ns,ns+1,ns+2);

    if(count == 6){
      // Change index start from 1 to 0
      n--;
      for(int i=0;i<3;i++) ns[i]--;
      if(n>=neighbours.size()){
        neighbours.resize(n+1);
        layout2d.resize(n+1);
      }
      neighbours[n] = vector<node_t>(ns,ns+3);
      layout2d[n] = coord2d(x,y);
    } else {    // Try
      char dummy[0x300];
      count = sscanf(line,"%d %s %s %d %d %d",&n,dummy,dummy,ns,ns+1,ns+2);
      if(count == 6){
        n--;
        for(int i=0;i<3;i++) ns[i]--;
        if(n>=neighbours.size())
          neighbours.resize(n+1);

      neighbours[n] = vector<node_t>(ns,ns+3);
      } else {
        //   fprintf(stderr,"Skipped line: %s\n",line);
      }
    }
  }
  N = neighbours.size();
}

// parse house of graphs
CubicGraph::CubicGraph(const unsigned int index, FILE* file){
  const int header_size = 15;

  // Get file size
  fseek(file, 0, SEEK_END);
  size_t file_size = ftell(file);

  //find number of vertices per graph
  //this only works for files with graphs of the equal size
  fseek(file, header_size, SEEK_SET);

  // Read the number N of vertices per graph.
  fread(reinterpret_cast<char*>(&N), 1, 1, file);
  if(N == 0){
    fread(reinterpret_cast<char*>(&N), 2, 1, file);
  }

  neighbours.resize(N);

  //only for files with graphs of the equal size
  unsigned int step;
  if(N<=255)
    {step = N * 4 + 1;}
  else
    {step = N * 8 + 3;}
  size_t address = header_size + step * index;

  //check if selected graphnumber is valid
  unsigned int graphs_per_file = (file_size - header_size ) /step;
  if(graphs_per_file -1 < index){
    cerr << "You asked for the " << index+1 << "th fullerene, but there are only " << graphs_per_file << " stored in this file." << std::endl;
    abort();
  }

  //the actual parsing of the selected graph:
  //go to the beginning of the selected graph
  fseek(file, address+1, SEEK_SET);

  if(N<=255){
    for(node_t u=0; u<N; ++u){
      for(int neighbour=0; neighbour<3; ++neighbour){
        unsigned char v;
        fread(reinterpret_cast<char*>(&v), 1, 1, file);
        neighbours[u].push_back(v-1);
      }
      //skip one byte because there is no interesting content
      fseek(file, 1, SEEK_CUR);
    }
  } else{
    fseek(file, 2, SEEK_CUR);//because three bytes are not read
    for(node_t u=0; u<N; ++u){
      for(int neighbour=0; neighbour<3; ++neighbour){
        unsigned short v;
        fread(reinterpret_cast<char*>(&v), 2, 1, file);
        neighbours[u].push_back(v-1);
      }
      fseek(file, 2, SEEK_CUR);
    }
  }

  N = neighbours.size();
}

// the jumps start at 0
// n is the number of vertices
CubicGraph::CubicGraph(const int n, const vector<int>& spiral_string, const jumplist_t& jumps) : CubicGraph() {
  assert(spiral_string.size() == n/2 + 2);

  Triangulation dual(spiral_string,jumps);
  CubicGraph G(dual.dual_graph());

  *this = G;
}

CubicGraph::CubicGraph(const full_spiral_name &fsn){
  assert(fsn.construction_scheme == full_spiral_name::CUBIC);
cout << "creating dual" << endl;
  Triangulation T(fsn);
cout << "dual created" << endl;
  *this = T.dual_graph();
}


bool CubicGraph::get_spiral_from_cg(const node_t f1, const node_t f2, const node_t f3, vector<int> &spiral, jumplist_t &jumps, const bool general) const
{
  spiral.clear();
  jumps.clear();
  vector<node_t> permutation_dummy;

  Triangulation Dual(this->dual_graph(6));

  if(!Dual.get_spiral(f1, f2, f3, spiral, jumps, permutation_dummy, general)) return false;
  assert(spiral.size()==N/2+2);
  return true;
}

bool CubicGraph::get_spiral_from_cg(vector<int> &spiral, jumplist_t &jumps, const bool canonical, const bool general, const bool pentagon_start) const
{
  spiral.clear();
  jumps.clear();

  Triangulation Dual(this->dual_graph(6));

  if(!Dual.get_spiral(spiral, jumps, general, pentagon_start)) return false;
  assert(spiral.size()==N/2+2);
  return true;
}


// Creates the (k,l)-Goldberg-Coxeter construction C_{(k^2+kl+l^2)n} of the current C_n
CubicGraph CubicGraph::GCtransform(const unsigned k, const unsigned l, const bool do_layout) const
{
  assert(layout2d.size()==N);
  Triangulation t(dual_graph());
  t.layout2d = t.tutte_layout(); // FIXME remove because unnecessary?
  Triangulation t_inflated(t.GCtransform(k,l));
  CubicGraph cg(t_inflated.dual_graph());
  return cg;
}

// returns the intersection of three sets, assuming that it's exactly one element (aimed at one vertex shared between three faces)
// yes, I know this is inefficient, but the vectors are short
int intersection_3(vector<int> &f1, vector<int> &f2, vector<int> &f3){
  for (int i=0; i<f1.size(); i++){
    for (int j=0; j<f2.size(); j++){
      if (f1[i] != f2[j]) continue;
      for (int k=0; k<f3.size(); k++){
        if (f1[i] == f3[k]) return f1[i];
      }
    }
  }
  assert(false);
  return -1;
}


// permutation of cub vertex numbers (ie, replace v by vertex_numbers[v], to get numbered vertices)
// T is the dual of the cubic graph
// perm is the permutation as returned by T.get_spiral()
// locants are vertices that should have small vertex numbers (as far as permitted by symmetry equivalent canonical spirals)
vector<int> CubicGraph::vertex_numbers(const Triangulation &T, const vector<vector<node_t>> &permutations, const vector<node_t> &locants) const{
  vector<tri_t> tri_faces = T.compute_faces();
  vector<node_t> vertex_numbers(N);
  vector<node_t> vertex_numbers_inv(N,INT_MAX);

  for(int p=0; p<permutations.size(); p++){
    const vector<node_t> &perm=permutations[p];
    // cout << "permutation of vertex numbers of triangulation: " << perm << endl;
    vector<int> perm_inv(perm.size());
    for(int i=0; i< perm.size(); i++) perm_inv[perm[i]] = i;
    // cout << "inverse permutation of vertex numbers of triangulation: " << perm_inv << endl;

    // cout << "triangles, relative to trig internal vertex numbers: " << tri_faces << endl;
    for(int i=0; i<tri_faces.size(); i++){
      for(int j=0; j<3; j++){
        tri_faces[i][j] = perm_inv[tri_faces[i][j]];
      }
    }
    // sort triples, first the vertices per face, then the faces
    for(tri_t& t: tri_faces){
      if(t[0] > t[1]) swap(t[0], t[1]);
      if(t[1] > t[2]) swap(t[1], t[2]);
      if(t[0] > t[1]) swap(t[0], t[1]);
    }
    sort(tri_faces.begin(), tri_faces.end());
    // cout <<  "triangles, relative to canon vertex numbers: " << tri_faces << endl;

    // permute back
    for(int i=0; i<tri_faces.size(); i++){
      for(int j=0; j<3; j++){
        tri_faces[i][j] = perm[tri_faces[i][j]];
      }
    }
    // cout << "triangles, in terms of internal numbers, sorted according to canon numbers: " << tri_faces << endl;

    // and then find the vertex in cub, which is part of faces {0,1,2},{0,1,5}, ...
    vector<face_t> cub_faces = compute_faces_oriented();
    // cout << "faces " << cub_faces << endl;
    vector<int> vertex_numbers_inv_tmp(N);
    for(int v=0; v<N; v++){
      face_t f0 = cub_faces[tri_faces[v][0]];
      face_t f1 = cub_faces[tri_faces[v][1]];
      face_t f2 = cub_faces[tri_faces[v][2]];
      vertex_numbers_inv_tmp[v] = intersection_3(f0, f1, f2);
    }
    // copy to vertex_numbers_inv?
    if(locants.size()==0){
      vertex_numbers_inv = vertex_numbers_inv_tmp;
      break;
    }
    for(int l=0; l<locants.size(); l++){
      if(vertex_numbers_inv_tmp[locants[l]] > vertex_numbers_inv[locants[l]]) break;
      if(vertex_numbers_inv_tmp[locants[l]] < vertex_numbers_inv[locants[l]]){
        vertex_numbers_inv = vertex_numbers_inv_tmp;
        break;
      }
    }
  }
  for(int i=0; i<vertex_numbers.size(); i++) vertex_numbers[vertex_numbers_inv[i]] = i;
  // cout << "permutation of cub vertex numbers (ie, replace v by vertex_numbers[v], to get numbered vertices): " << vertex_numbers << endl;
  return vertex_numbers;
}

