#include "spherical-harmonic.hh"
#include "y3table.hh"

namespace RealSphericalHarmonic {

  double Y3(int l, int m, const coord3d& u)
  {
    size_t index  = l*l + l + m;
    size_t length = Y3Table::indices[index+1]-Y3Table::indices[index];

    double sum = 0;
    for(int i=0;i<length;i++){
      const uint16_t       &key(Y3Table::keys  [index+i]);
      const double &coefficient(Y3Table::values[index+i]);
      
      int n0 = key&0x1f, n1 = (key>>5)&0x1f, n2 = (key>>10)&0x1f;
      sum += coefficient * pow(u[0],n0) * pow(u[1],n1) * pow(u[2],n2);
    }
    return sum;
  }

  vector<Ylm_coefficient> decompose_polyhedron(const Polyhedron& P, int Lmax)
  {
    const vector<coord3d> &xs(P.points);
    int Nx = P.points.size(), Nlm = (Lmax+1)*(Lmax+1);
    
    vector<double> norm(Nlm);
    vector<double> sum(Nlm);
    for(int i=0;i<Nx;i++){
      const double  r = xs[i].norm();
      const coord3d u = xs[i]/r;

      for(int l=0;l<=Lmax;l++)
        for(int m=-l;m<=l;m++){
          const double Ylm = Y3(l,m,u);
          norm[l*l+l+m] += Ylm*Ylm;
          sum[l*l+l+m] += Ylm*r;
        }
    }
    
    vector<Ylm_coefficient> C(Nlm);
    for(int l=0,i=0;l<=Lmax;l++)
      for(int m=-l;m<=l;m++,i++){
        C[i].l = l;
        C[i].m = m;
        C[i].coefficient = sum[i]/sqrt(norm[i]);
      }
    return C;
  }
  
}
