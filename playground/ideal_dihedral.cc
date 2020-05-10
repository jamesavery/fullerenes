#include <stdio.h>
#include <stdlib.h>
#include "fullerenes/geometry.hh"


int main(int ac, char **av)
{
  vector<double>   l(2*2);	// Edge lengths. l[i][j] corresponds to edge-length between two faces of size (5+i) and (5+j).
  vector<double> dih(2*2*2); 	// Ideal dihedrals dih[i][j][k] corresponds to face sizes (5+i), (5+j), and (5+k).

  if(ac<=3){
    fprintf(stderr,"Syntax: %s l_pp l_hp l_hh \n\
(The lengths for hexagon-hexagon edges, hexagon-pentagon, and pentagon-pentagon)\n",
	    av[0]);
    return -1;
  }

  double l_pp = strtod(av[1],0);
  double l_hp = strtod(av[2],0);
  double l_hh = strtod(av[3],0);

  l = {l_pp,l_hp,
       l_hp,l_hh};

  cout << "lengths = np.array(" << l << ").reshape(2,2)\n";
  
//        t   B   s      
//          \   /
//        C   u   A
//            |
//            r  
  for(int i=0;i<2;i++)
    for(int j=0;j<2;j++)
      for(int k=0;k<2;k++)
	dih[i*2*2+j*2+k] = coord3d::ideal_dihedral(5+i,5+j,5+k,
						   l[i*2+k], l[j*2+i], l[k*2+j]);

  cout << "dihedrals = np.array("<<dih << ").reshape(2,2,2)\n";
  
  return 0;
}
