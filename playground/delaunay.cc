#include "libgraph/fullerenegraph.hh"
#include "libgraph/triangulation.hh"
#include <math.h>

using namespace std;

void insert_before(vector<int> &v, const int before_what, const int value){
  vector<int>::iterator pos = std::find(v.begin(), v.end(), before_what);
  v.insert(pos, value);
}

Triangulation Delaunayify(Triangulation T, double distances[12][12]){

//         C                     C       
//     c / | \ b               /   \      .
//     /   |   \             /       \    .
//   D     |     B   --->  D-----------B 
//     \   |   /             \       /   
//     d \ | / a               \   /     
//         A                     A       
 

  int A, B, C, D;
  unsigned int delaunay_flips = 0;

  auto flip = [&](){
    delaunay_flips++;
    T[A].erase(find(T[A].begin(), T[A].end(), C));
    T[C].erase(find(T[C].begin(), T[C].end(), A));
    insert_before(T[B], A, D);
    insert_before(T[D], C, B);
  };

  while(delaunay_flips != 0){
    delaunay_flips=0;
    for(node_t u=0; u<T.N; ++u){
      for(int j=0; j<T[u].size(); ++j){

        A = u;
        C = T[u][j];
        B = T[u][(j+1)%T[u].size()];
        D = T[u][(j+2)%T[u].size()];
        const double a = distances[A][B],
                     b = distances[B][C],
                     c = distances[C][D],
                     d = distances[D][A],
                     AC = distances[A][C];
        const double beta = acos((a*a + b*b - AC*AC) / 2*a*b),
                     delta = acos((c*c + d*d - AC*AC) / 2*c*d);
        if(beta + delta > 180.0) flip();
      }
    }
    cout << "flips in this sweep: " << delaunay_flips << endl;
  }

  return T;
}


int main(){

cout << "foo" << endl;

return 0;
}



