#version 3.6;
#include "textures.inc"
#include "colors.inc"

background { White }

camera {
  location <0,0,-15>
  look_at <0,0,0>
  focal_point <0,0,0>
  up <0,1,0>
  right <1,0,0>
}

light_source {
  <-15,-15,-15>
  color rgb <1,1,1>
}

light_source {
  <15,15,15>
  color rgb <1,1,1>
}

light_source {
  <-15,15,-15>
  color rgb <1,1,1>
}

union {
#declare i=0;
#while (i<dimension_size(tris,1))
  #declare k=tris[i][0];
  #declare l=tris[i][1];
  #declare m=tris[i][2];
  polygon {
    4,
    <cpoints[k][0],cpoints[k][1],cpoints[k][2]>,
    <cpoints[l][0],cpoints[l][1],cpoints[l][2]>,
    <cpoints[m][0],cpoints[m][1],cpoints[m][2]>,
    <cpoints[k][0],cpoints[k][1],cpoints[k][2]>

    pigment { facecolour transmit (1.0-faceopacity) }
    finish { phong .8 reflection {0.3} }
  }
  #declare i=i+1;
#end
}

union {
#declare i=0;
#while (i<Nvertices)
  sphere {
    <layout3D[i][0],layout3D[i][1],layout3D[i][2]>, nodediameter
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
    <layout3D[uu][0],layout3D[uu][1],layout3D[uu][2]>,
    <layout3D[vv][0],layout3D[vv][1],layout3D[vv][2]>,
    edgewidth

    texture {
      pigment { edgecolour }
      finish { phong .5 }
    }
  }
  #declare i=i+1;
#end
}
