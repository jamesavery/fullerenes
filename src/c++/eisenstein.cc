#include "fullerenes/eisenstein.hh"

// (-1,1)  x___x (0,1)
//        / \ / \ 
// (-1,0)x---x---x (1,0)
//        \ / \ /
//  (0,-1) x---x (1,-1)
Eisenstein Eisenstein::unit[7] = {{1,0},{0,1},{-1,1},{-1,0},{0,-1},{1,-1},{1,0}};

				  
