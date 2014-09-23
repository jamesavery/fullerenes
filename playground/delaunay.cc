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
    T.neighbours[A].erase(find(T.neighbours[A].begin(), T.neighbours[A].end(), C));
    T.neighbours[C].erase(find(T.neighbours[C].begin(), T.neighbours[C].end(), A));
    insert_before(T.neighbours[B], A, D);
    insert_before(T.neighbours[D], C, B);
  };

  while(delaunay_flips != 0){
    delaunay_flips=0;
    for(int i=0; i<T.neighbours.size(); ++i){
      for(int j=0; j<T.neighbours[i].size(); ++j){

        A = i;
        C = T.neighbours[i][j];
        B = T.neighbours[i][(j+1)%T.neighbours[i].size()];
        D = T.neighbours[i][(j+2)%T.neighbours[i].size()];
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



