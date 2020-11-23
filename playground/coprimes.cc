#include <stdio.h>
#include <iostream>
#include <chrono>
#include <vector>
#include <cmath>

using namespace std;
typedef uint32_t Int;

vector<pair<Int,Int>> Farey_sequence(Int Nmax)
{
  Int length_bound = 1+Nmax+int(std::floor(Nmax*Nmax)*0.305);
  vector<pair<Int,Int>> result;  
  result.reserve(length_bound);

  Int a = 0, b = 1, c = 1, d = Nmax;
  result.push_back({a,b});
  while(c <= Nmax){
    Int k = (Nmax+b) / d;
    //    printf("(a,b,c,d) = (%d,%d,%d,%d), k = %d\n",a,b,c,d,k);    
    Int a_next = c;
    Int b_next = d;
    c = k*c - a;
    d = k*d - b;
    a = a_next;
    b = b_next;
    result.push_back({a,b});
  }
  return result;
}


void Farey_inplace(Int Nmax, vector<pair<Int,Int>> &result)
{
  result.clear();
  
  Int a = 0, b = 1, c = 1, d = Nmax;
  result.push_back({a,b});
  while(c <= Nmax){
    Int k = (Nmax+b) / d;
    //    printf("(a,b,c,d) = (%d,%d,%d,%d), k = %d\n",a,b,c,d,k);    
    Int a_next = c;
    Int b_next = d;
    c = k*c - a;
    d = k*d - b;
    a = a_next;
    b = b_next;
    result.push_back({a,b});
  }
}


pair<int,int> Farey(Int Nmax)
{
  int a = 0, b = 1, c = 1, d = Nmax;

  while(c <= Nmax){
    Int k = (Nmax+b) / d;
    //    printf("(a,b,c,d) = (%d,%d,%d,%d), k = %d\n",a,b,c,d,k);    
    Int a_next = c;
    Int b_next = d;
    c = k*c - a;
    d = k*d - b;
    a = a_next;
    b = b_next;
  }
  return {a,b};
}

int main()
{
  typedef std::chrono::high_resolution_clock Clock;
  Int Nmax = 1000;
  double length_bound = std::floor(1+Nmax+Nmax*Nmax*0.305);
  cout << Nmax << ", " << int(length_bound) << endl;
  
  vector<pair<Int,Int>> result;  
  result.reserve(int(length_bound));
  

  auto t1 = Clock::now();
  for(int i=0;i<1000;i++)
    Farey(Nmax);
  auto t2 = Clock::now();
  auto ns = chrono::duration_cast<chrono::nanoseconds>(t2-t1);  
  cout << (ns.count()/(1000*length_bound)) << " nanoseconds.\n";
}
