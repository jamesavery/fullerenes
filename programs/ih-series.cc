// 5(6)^a (5(6)^b)^4 5(6)^c (5(6)^d)^4 5(6)^e 5
#include <tuple>
#include <array>
#include <vector>
#include <iostream>
#include <assert.h>

using namespace std;

// template <typename T> ostream& operator<<(ostream& s, const vector<T>& v) 
// { 
//   for(int i=0;i<v.size();i++) s<< v[i] << (i+1<v.size()? ",":"");
//   return s; 
// }


template <typename T> ostream& operator<<(ostream& s, const vector<T>& v) 
{ 
  for(int i=0;i<v.size();i++) s<< v[i] << (i+1<v.size()? ",":"");
  return s; 
}

array<int,5> IhGC_to_abcde(int i, int j)
{
  int a,b,c,d,e;
  a = ((i+j)*(i+j)*5 - 5*i - 3*j)/2 - 1 + i*(j==0); // Eq. (3)
  b = i+j-1;					    // Eq. (4)
  c = (5*i+1)*(i-1) + j*(5*i-3);		    // Eq. (5)
  d = b;
  e = (5*i*i + 15*j*j - 3*i - 7*j)/2 - i*(j==0);    // Eq. (6)
 
  assert((a + 4*b + c + 4*d + e) == 10*(i*i+i*j+j*j-1)); // Eq. (7)
  
  return {a,b,c,d,e};
}

vector<int> abcde_to_rspi(array<int,5> abcde)
{
  vector<int> rspi(12);
  // Decoding Equation (2)
  auto [a,b,c,d,e] = abcde;
    
  int Nf = a+4*b+c+4*d+e+12;

  rspi[0]  = 1; 		 // 5 
  rspi[1]  = rspi[0]+(a+1);	 // 6^a 5
  rspi[2]  = rspi[1]+(b+1);	 // 6^b 5
  rspi[3]  = rspi[2]+(b+1);	 // 6^b 5 
  rspi[4]  = rspi[3]+(b+1);      // 6^b 5 
  rspi[5]  = rspi[4]+(b+1);      // 6^b 5
  rspi[6]  = rspi[5]+(c+1);      // 6^c 5
  rspi[7]  = rspi[6]+(d+1);      // 6^d 5
  rspi[8]  = rspi[7]+(d+1);      // 6^d 5
  rspi[9]  = rspi[8]+(d+1);      // 6^d 5
  rspi[10] = rspi[9]+(d+1);      // 6^d 5
  rspi[11] = rspi[10]+(e+1);      // 6^e 5

  if(!(rspi[11] == Nf)){
    cout << "rspi[11] = " << rspi[11] << endl;
  }

  return rspi;
}

int main(int ac, char **av)
{
  assert(ac>=2);
  
  int I = strtol(av[1],0,0);
  int J = strtol(av[2],0,0);


  for(int i=1;i<=I;i++)
    for(int j=0;j<=i && j<=J;j++){
      array<int,5> abcde = IhGC_to_abcde(i,j);
      vector<int> rspi  = abcde_to_rspi(abcde);
      auto [a,b,c,d,e] = abcde;
      int Nf = a+4*b+c+4*d+e+12;
      int N  = (Nf-2)*2;
      cout << "C"<<N<<"-[GS:"<<rspi<<"]-fullerene\n";
    }      
  
  return 0;
}
