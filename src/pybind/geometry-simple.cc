#include <chrono>
#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <math.h>
using namespace std;

#include <array>
#include <tuple>
typedef double real_t;
typedef array<real_t,3> coord3d;

constexpr int N = 400;
typedef array<coord3d, N> VectorNx3;
typedef array<int, N*3>   CubicArcs;


__attribute__((always_inline))  inline real_t dot(const coord3d& x, const coord3d& y)
{
  return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

__attribute__((always_inline))  inline pair<real_t,coord3d> split_norm(const coord3d& x)
{
  real_t r = sqrt(dot(x,x));

  return {r, {x[0]/r, x[1]/r, x[2]/r}};
}



int main(int ac, char **av)
{
  array<real_t,3> x[1000];

  for(int i=0;i<1000;i++) x[i] = {double(i),double(i),double(i)};

  real_t r;
  coord3d xhat;  
  real_t  r_sum = 0;

  for(int i=0;i<1000;i++){
    tie(r,xhat) = split_norm(x[i]);
    r_sum += r;
  }

  return int(r_sum);
}
