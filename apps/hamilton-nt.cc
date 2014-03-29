// Calculates number of Hamilton cycles in a D5{h,d}-nanotube.
#include <stdio.h>
#include <vector>
#include <string>

using namespace std;

class D5 {
public:
  typedef enum {A,B} start_type;
  
  vector<int> shaded_data;
  
  bool  shaded(int i,int j) const { return shaded_data[i] & (1<<(j%5)); }
  int&  shaded(int i)             { return shaded_data[i]; }

  bool safe(int i) const {
    // Ring has shaded face adjacent to every vertex
    bool reaches_all = true;
    // No shaded face has both "parents" shaded.
    bool occludes_vertex = false;

    for(int j=0;j<5;j++){
      reaches_all     &= shaded(i,j) || shaded(i,j+1) || shaded(i-1,j);
      occludes_vertex |= shaded(i,j)   && (shaded(i-1,j) && shaded(i-1,j+1));
      occludes_vertex |= shaded(i-1,j) && (shaded(i,j) && shaded(i,j+1));
    }
    return reaches_all && !occludes_vertex;
  }


  D5(int n, start_type start) : shaded_data(n) {
    switch(start){
    case A:
      shaded(0) = 22;//10110
      break;
    case B:
      shaded(0) = 15;//01111
      break;
    }
  }

  static string binary(int w) {
    char result[5];
    for(int i=0;i<5;i++) result[5-i-1] = (w>>i)&1? '1':'0';
    return string(result,5);
  }


  bool is_hamiltonian() const 
  {
    return true;		// hrow to test?
  }

  size_t count_hamilton_paths(int i=1) 
  {
    for(int k=0;k<i;k++) printf(" ");
    printf("ring %d = %s\n",i-1,binary(shaded(i-1)).c_str());
    if(i == shaded_data.size()) return is_hamiltonian();

    size_t cnt = 0;
    for(shaded(i)=0;shaded(i)<32;shaded(i)++)
      if(safe(i)) cnt += count_hamilton_paths(i+1);
    return cnt;
  }
};


int main()
{
  D5 nt(4,D5::A);

  printf("%ld\n",nt.count_hamilton_paths());
}
