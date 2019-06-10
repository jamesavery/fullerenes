#version 3.6;
#include "textures.inc"
#include "colors.inc"

background { White }

camera {
  location <0,0,-11>
  look_at <0,0,0>
  focal_point <0,0,0>
  up <0,1,0>
  right <1,0,0>
}

light_source {
  <-10,-10,-10>
  color rgb <1,1,1>
}

light_source {
  <-10,-10,10>
  color rgb <1,1,1>
}


#declare i=0;
#while (i<Nvertices)
  sphere {
    <layout2D[i][0],layout2D[i][1],0>, nodediameter
    texture {
      pigment { nodecolour }
      finish { phong 1 }
    }
  }
  #declare i = i+1;
#end

#declare i=0;
#while (i<Nedges)
  #declare uu=edges[i][0];
  #declare vv=edges[i][1];
  cylinder {
    <layout2D[uu][0],layout2D[uu][1],0>,
    <layout2D[vv][0],layout2D[vv][1],0>, edgewidth

    texture {
      pigment { edgecolour }
      finish { phong .5 }
    }
  }
  #declare i=i+1;
#end
