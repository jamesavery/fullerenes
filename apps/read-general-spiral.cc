#include <string>

#include "libgraph/triangulation.hh"
#include "libgraph/spiral.hh"


map<string,string> spiral_paper_examples{{
    {"Tetrahedron", "[3,3,3,3]"},
    {"Truncated_tetrahedron","[CS: 6,3,6,3,6,3,6,3]"},
    {"Truncated_tetrahedron GS","[6,1; 3,6,6,6,3,3,6,3]"},
    {"Triakis_tetrahedron","[T,CS: 6,3,6,3,6,3,6,3]"},
    {"Omnitruncated_octahedron_spiral","[4, 6,6,6, 6,4,6,4,6,4,6,4,6,4]"},
    {"Omnitruncated_octahedron_full",  "Oh-[4, 6,6,6, 6,4,6,4,6,4,6,4,6,4]-24-cage"},
    {"Tutte_graph","[CS:11, 1, 17, 1; 5, 10, 5, 5, 5, 9, 5, 4, 5, 4, 4, 5, 4, 10, 5, 5, 5, 5, 5, 10, 4, 5, 5, 4, 5]"},
    {"Tutte_molecule","C3–[CS: 11,1,17,1; 5,10,5,5,5,9,5,4,5,4,4,5,4,10,5,5,5,5,5,10,4,5,5,4,5]-C46-cage"},
      // Different ways of writing the Td-C100
    {"Td-C100_shortest","[2,8,9,23,24,28,29,37,41,45,46,52]-fullerene"},
    {"Td-C100_full_GS","Td-[GS: 43,2; 1,4,5,26,27,31,32,40,43,47,48,52]-C100-fullerene"},
    {"Td-C100_roid_a","[GS: 43,2; 1,4,5,26,27,31,32,40,43,47,48,52]-(5)-C100-fulleroid"},
    {"Td-C100_roid_b","[GS: 43,2; 1,4,5,26,27,31,32,40,43,47,48,52]-(5)6-fulleroid"},
    {"Td-C100_roid_c","[GS: 43,2; 1,4,5,26,27,31,32,40,43,47,48,52]-(5)_6-fulleroid"},
    {"Td-C100_roid_d","[GS: 43,2; ; 1,4,5,26,27,31,32,40,43,47,48,52]-(4,5)-fulleroid"},
    {"Td-C100_roid_e","[GS: 43,2; ; 1,4,5,26,27,31,32,40,43,47,48,52; ]-(4,5,8)-fulleroid"},
      // Non-spiralable fullerene examples
    {"NS-T-C380","[GS: 162,2,186,3; 1,4,5,126,127,136,137,155,162,171,172,186]-C380-fullerene"},
    {"NS-D3-C384","D3-[GS:171,8,178,9; 1,3,4,5,168,169,170,178,190,191,192,194]-C384-fullerene"},
    {"NS-D3-C440","D3-[CS:62,1;39,40,41,62,170,171,197,198,218,219,220,222]-C440-fullerene"},
    {"NS-D3-C672","D3-[GS:297,8,310,9; 1,10,12,14,260,262,264,310,324,326,328,338]-C672-fullerene"},
      // Non-cubic polyhedra nanoparticle-examples
    {"M12L24","Oh-[LF,GS: 9,1,12,1,15,1,24,3; "
	   "4, 8,8,8,8, "
	   "3, 8, 4, 3, 8, 4, 3, 8, 4, 3, 8, 4,"
 	   "8, 3,8, 3,8, 3,"
  	   "3, 8, 4]-12-cage"},
    {"M24L48","Oh-[LF,GS: 22,1,27,1,42,3; "
	   "3, 8,8,8, "
	   "4,8, 4,8, 4,8, 4,8, 4,8, 4,8, "
	   "3,8,4,8,4, 3,8,4,8,4, 3,8,4,8,4, "
	   "8,4,3, "
	   "8,4, 8,4, 8,4, 8,4, "
	   "8,3"
	"]-24-cage"},
    {"M30L60","O-[LF,GS: 54,1,57,1,59,1; "
	   "3,8,8,8, "
	   "4,8, 4,8, 4,8, 4,8, 4,8, 4,8, 4,8, "
	   "3,8,4,8,4,8, 3,8,4,8,4,8, 3,8,4,8,4,8, "
	   "4,8, "
	   "3,8,4,8,4,8, 3,8,4,8,4,8, "
	   "8,3, "
	   "8,4,4, 8,4,4, 8,4,4, "
	   "8,3"
	"]-30-cage"}
  }};


int main(int ac, char **av)
{
  for(auto example: spiral_paper_examples){
    cout << example.first << "_name = \"" << example.second << "\";\n";
    cout << example.first << " = " << full_spiral_name(example.second) << ";\n\n";
  }

  return 0;
}
