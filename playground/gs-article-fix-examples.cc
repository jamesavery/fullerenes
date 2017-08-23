void fix_up_C1C100_layout(PlanarGraph& g)
{
  double radius0 = 1, radius1 = 0.8;
  face_t outer_face{{55,99,65,89,91}};
  outer_face = outer_face+(-1);
  
  vector<coord2d> &xs(g.layout2d);

  for(int i=0;i<g.N;i++)
    xs[i] = (xs[i]+coord2d{0,-.02})*2.5;
      
  for(int i=0;i<outer_face.size();i++)
    xs[outer_face[i]] =
      coord2d{sin(2*M_PI*i/outer_face.size()),cos(2*M_PI*i/outer_face.size())}*radius0;

  // Second-to-outer face face
  xs[90-1] = xs[55-1]*.9;
  xs[26-1] = xs[99-1]*.9;
  xs[98-1] = xs[65-1]*.9;
  xs[92-1] = xs[89-1]*.9;
  xs[28-1] = xs[91-1]*.9;

  xs[1-1]  = (xs[90-1]+xs[26-1])/2.0;
  
  xs[4-1]  = (xs[90-1]*2+xs[28-1]  )/3;
  xs[23-1] = (  xs[90-1]+xs[28-1]*2)/3;

  xs[84-1]  = (xs[28-1]*2+  xs[92-1])/3;
  xs[82-1] = (   xs[28-1]+xs[92-1]*2)/3;

  xs[50-1]  = (xs[92-1]*2+xs[98-1]  )/3;
  xs[88-1] = (   xs[92-1]+xs[98-1]*2)/3;

  xs[96-1]  = (xs[98-1]*2+xs[26-1]  )/3;
  xs[24-1] = (   xs[98-1]+xs[26-1]*2)/3;

  // Third-outer ring
  xs[54-1] = (xs[4-1]+xs[1-1])/2;
  xs[ 2-1] = (xs[1-1]+xs[24-1])/2;  
  xs[93-1] = xs[23-1]*.9;
  xs[22-1] = xs[84-1]*.9;
  xs[30-1] = xs[82-1]*.9;
  xs[79-1] = xs[50-1]*.9;
  xs[83-1] = xs[96-1]*.9;
  xs[61-1] = xs[ 1-1]*.9;

  xs[67-1] = xs[92-1]*.95;

  xs[33-1] = (xs[54-1]*5+xs[93-1]  )/6;
  xs[97-1] = (xs[54-1]*4+xs[93-1]*2)/6;

  xs[100-1] = (xs[ 2-1]*5+xs[83-1]  )/6;
  xs[ 14-1] = (xs[ 2-1]*4+xs[83-1]*2)/6;
  
  xs[ 94-1] = (xs[97-1]+xs[14-1])/2;

  xs[37-1] = (xs[37-1]+xs[94-1]*5)/5;
  xs[ 4-1] = (xs[ 4-1]+xs[90-1]*2)/3;
  xs[24-1] = (xs[24-1]+xs[26-1]*2)/3;  

  xs[ 6-1] = (xs[93-1]*2+xs[37-1]  )/3;
  xs[73-1] = (xs[93-1]  +xs[37-1]*2)/3;

  xs[21-1] = (xs[83-1]*2+xs[37-1]  )/3;
  xs[45-1] = (xs[83-1]  +xs[37-1]*2)/3;

  xs[54-1] = (xs[33-1]  +xs[ 4-1]*2)/3;
  xs[ 2-1] = (xs[100-1] +xs[24-1]*2)/3;    

  xs[60-1] = (xs[60-1]+xs[45-1]*5)/5;
  xs[76-1] = (xs[76-1]+xs[73-1]*5)/5;
  xs[22-1] = (xs[84-1]+xs[6-1]*5)/5;
  xs[86-1] = (xs[88-1]+xs[21-1]*5)/5;

  xs[ 6-1] = (xs[93-1]+xs[22-1])/2;
  xs[21-1] = (xs[86-1]+xs[83-1])/2;

  xs[68-1] = (xs[64-1]+xs[67-1]*5)/6;
  xs[84-1] = (xs[84-1]+xs[28-1]*3)/4;
  xs[88-1] = (xs[88-1]+xs[98-1]*3)/4;    
  xs[82-1] = (xs[82-1]+xs[84-1]*3)/4;
  xs[50-1] = (xs[50-1]+xs[88-1]*3)/4;

  xs[22-1] = (xs[ 6-1]+xs[84-1]  )/2;  
  xs[30-1] = (xs[82-1]*3+xs[22-1]  )/4;
  xs[70-1] = (xs[82-1]*2+xs[22-1]*2)/4;
  xs[74-1] = (xs[82-1]  +xs[22-1]*3)/4;

  xs[86-1] = (xs[21-1]+xs[88-1]  )/2;  
  xs[78-1] = (xs[86-1]*3+xs[50-1]  )/4;
  xs[77-1] = (xs[86-1]*2+xs[50-1]*2)/4;
  xs[79-1] = (xs[86-1]  +xs[50-1]*3)/4;
  
  xs[63-1] = (xs[63-1]  +xs[77-1]*5)/6;
  xs[69-1] = (xs[69-1]  +xs[70-1]*5)/6;  
  xs[11-1] = (xs[74-1]*2+xs[ 6-1]+xs[76-1]*2)/5;
  xs[44-1] = (xs[78-1]*2+xs[21-1]+xs[60-1]*2)/5;  

  xs[66-1] = (xs[69-1]  +xs[68-1]*7)/8;
  xs[19-1] = (xs[63-1]  +xs[68-1]*7)/8;  

  xs[ 7-1] = (xs[ 7-1]  +xs[66-1]*2)/3;
  xs[81-1] = (xs[81-1]  +xs[19-1]*2)/3;
  xs[17-1] = (xs[11-1]*2+xs[69-1]  )/3;
  xs[13-1] = (xs[11-1]  +xs[69-1]*2)/3;
  xs[27-1] = (xs[44-1]*2+xs[63-1]  )/3;
  xs[36-1] = (xs[44-1]  +xs[63-1]*2)/3;

  xs[48-1] = (xs[ 7-1]*3+xs[13-1]  )/4;
  xs[ 5-1] = (xs[ 7-1]  +xs[13-1]*3)/4;
  xs[ 8-1] = (xs[81-1]*3+xs[36-1]  )/4;
  xs[53-1] = (xs[81-1]  +xs[36-1]*3)/4;        

  xs[64-1] = (xs[76-1]+xs[60-1])/2;
  xs[ 3-1] = (xs[ 3-1]+xs[64-1]*8)/9;
  xs[15-1] = (xs[ 3-1]+xs[17-1]*4)/5;
  xs[85-1] = (xs[ 3-1]+xs[27-1]*4)/5;

  xs[ 9-1]  = (xs[15-1]*5+xs[85-1]  )/6;
  xs[59-1]  = (xs[15-1]*3+xs[85-1]*3)/6;
  xs[51-1]  = (xs[15-1]  +xs[85-1]*5)/6;    

  xs[40-1] = (xs[ 9-1]+xs[ 5-1]*2)/3;
  xs[32-1] = (xs[51-1]+xs[53-1]*2)/3;  

  xs[10-1] = (xs[40-1]*3+xs[59-1]  )/4;
  xs[57-1] = (xs[40-1]*2+xs[59-1]*2)/4;
  xs[75-1] = (xs[71-1]  +xs[59-1]*7)/8;  
  
  xs[62-1] = (xs[32-1]*3+xs[59-1]  )/4;
  xs[49-1] = (xs[32-1]*2+xs[59-1]*2)/4;

  xs[80-1] = xs[68-1]*.9;
  xs[25-1] = xs[80-1]*.9;  
  
  face_t fixed_vertices{{55,91,89,65,99,
	26,1,90,4,23,28,84,82,92,50,88,98,96,24,
	54,93,22,30,79,83,2,61,
	67,33,
	97,100,
	14,
	94,37,6,73,21,45,
	60,76,22,86,68,74,30,70,77,78,
	63,69,11,44,66,19,
	7,81,17,13,27,36,48,5,8,53,64,3,15,85,
	9,59,51,40,32,10,62,57,75,49,
	80,25}};
  fixed_vertices = fixed_vertices+(-1);
  xs = g.tutte_layout_iterative(fixed_vertices,xs);
}  

void fix_up_TdC100_layout(PlanarGraph& g)
{
  face_t outer_face{{72, 97, 99, 82, 89}};
  double radius0 = 0.5, radius1 = 0.43;

  for(int i=0;i<g.N;i++)
    g.layout2d[i] = (g.layout2d[i]+coord2d{0,.07})*1.6;
      
  for(int i=0;i<outer_face.size();i++)
    g.layout2d[outer_face[i]] =
      coord2d{sin(2*M_PI*i/outer_face.size()),cos(2*M_PI*i/outer_face.size())}*radius0;

  face_t next_ring{{84,56,88,57,32,86,85,39,23,63,0,42,1}};
  vector<int> nodes{{0,2,5,8,11}};
  vector<vector<int>> intermediates{{1},{3,4},{6,7},{9,10},{12}};
      
  for(int i=0;i<nodes.size();i++)
    g.layout2d[next_ring[nodes[i]]] = coord2d{sin(-2*M_PI*i/nodes.size()),
					      cos(2*M_PI*i/nodes.size())}*radius1;

  g.layout2d[56] = (g.layout2d[84]+g.layout2d[88])/2.0;

  g.layout2d[57] = (g.layout2d[88]*2+g.layout2d[86])/3.0;
  g.layout2d[32] = (g.layout2d[88]  +g.layout2d[86]*2)/3.0;

  g.layout2d[85] = (g.layout2d[86]*2+g.layout2d[23])/3.0;
  g.layout2d[39] = (g.layout2d[86]  +g.layout2d[23]*2)/3.0;
      
  g.layout2d[63] = (g.layout2d[23]*2+g.layout2d[42])/3.0;
  g.layout2d[0]  = (g.layout2d[23]  +g.layout2d[42]*2)/3.0;

  g.layout2d[1]  = (g.layout2d[42]+g.layout2d[84])/2.0;

  // Next-next layer
  g.layout2d[47] = g.layout2d[56]*0.85;
  g.layout2d[2] = g.layout2d[1]*0.85;
  g.layout2d[91] = (g.layout2d[2]+g.layout2d[47])/2.0;
}  
