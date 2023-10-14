

template <typename t> t *shalloc(uint8_t *smem, size_t nelem, size_t &offset)
{
  if(offset%sizeof(t)>0) offset += sizeof(t)-(offset%sizeof(t));
  t *new_ptr = reinterpret_cast<t*>(smem+offset);
  offset += nelem*sizeof(t);

  return t;
}

KERNEL_SIGNATURE(T,U)
void tutte_layout_(IsomerBatch<T,U> B, const size_t iterations)
{
  template typename T::real_t          real_t;
  template typename T::node_t          node_t;
  template typename coord2d<real_t>    coord;

  const size_t N = B.n_atoms;

  LOCAL_DEFINE(relative_change,real_t,N); // real_t relative_change_array[N];
  
  // HEADER_BOILERPLATE;
  // TODO: IsomerBatch -> potentiel een isomer; isomer iteration til boiler.
  for (int isomer_idx = blockIdx.x; isomer_idx < B.isomer_capacity; isomer_idx+= gridDim.x)
    if (B.statuses[isomer_idx] != IsomerStatus::EMPTY){
      node_t  outer_face[6];	
      node_t  outer_face_vertex   = 0;      
      uint8_t outer_face_len = G.get_face_oriented(0,G[0][0], outer_face);
      
      size_t isomer_offset = isomer_idx * N;
      CubicGraph G(&B.cubic_neighbours[isomer_offset*3]);

      size_t  shared_offset = 0;
      real_t *reduction_smem = shalloc<real_t>(sharedmem,Block_Size_Pow_2,shared_offset);
      coord  *xys            = shalloc<coord> (sharedmem,N,shared_offset);
      coord  *newxys         = shalloc<coord> (sharedmem,N,shared_offset);
      bool   *fixed_array    = shalloc<bool*> (sharedmem,outer_face_len, shared_offset);

      LOOP(i,N){ fixed_array[i] = false; }
      LOOP(i,outer_face_len){
	size_t fi = outer_face[i];	
	xys[fi] = {sin(a*2.0L*M_PI/real_t(outer_face_len)),
	           cos(a*2.0L*M_PI/real_t(outer_face_len))};
	newxys[fi]      = xys[fi];
	fixed_array[fi] =  true;	
      }
      
      bool converged    = false;

      for (size_t i = 0; i < iterations && !converged; i++){      
	real_t max_change = 0.0;
	
	LOOP(a,N){
	  LOCAL_USE(a,relative_change); // auto &relative_change = relative_change_array[a];
	  
	  node_t ns[3] = {G[a][0], G[a][1], G[a][2]};
	  xys[a]       = {0.0, 0.0};
	  
	  bool fixed            = fixed_array[a];
	  coord2 neighbour_mean = {0.0,0.0};    
	  real_t neighbour_dist = 0.0;
	  
	  for (uint8_t j = 0; j < 3; j++) neighbour_mean += xys[ns[j]];
	  neighbour_mean /= real_t(3.0L);
	  
	  // Calculate the new position of the point
	  if(!fixed) newxys[a] = xys[a]*real_t(0.15) + neighbour_mean*real_t(0.85L);

	  // Calculate the distance between neighbours
	  for (uint8_t j = 0; j < 3; j++)
	    neighbour_dist += norm(xys[a]-xys[ns[j]])/real_t(3.0L);
        
	  BLOCK_SYNC; // TODO: Skal denne bruges?
	  real_t relative_change = 0.0;
	  // Calculate the relative change
	  if (neighbour_dist > real_t(0.0) && !fixed)
            relative_change = norm(xys[a] - newxys[a])/neighbour_dist;
	}

	// Reduce the relative change to find the maximum change
	// TODO: Can we just drop the reduction, and iterate a fixed number of steps?
	BLOCK_SYNC;
	real_t max_change = reduction_max(sharedmem, relative_change);
	converged = max_change <= 100*numeric_limits<real_t>::epsilon();

	// Update the position of the point
	xys[a] = newxys[a];
      }
    }

  
  LOOP(a,N){
    output_xys[a] = xys[a];
  }
}

