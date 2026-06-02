#include <sstream>
#include "fullerenes/graph.hh"
#include "fullerenes/cubicgraph.hh"

struct ToleranceLess {
  const double tolerance;
  ToleranceLess(const double tolerance) : tolerance(tolerance) {}
  bool operator()(const coord2d& x,const coord2d& y) const { return x<y && (y-x).norm() > tolerance; }
};

vector<coord2d> PlanarGraphView::tutte_layout(node_t s, node_t t, node_t r, unsigned int face_max) const
{
  if(count_edges() == 0)       // empty graph
    return vector<coord2d>(); // -> empty layout

  if(s<0) s = 0;
  if(t<0) t = (*this)[s][0];

  face_t outer_face;
  if(is_consistently_oriented()){
    outer_face = get_face_oriented({s,t},face_max);
  } else {
    // For unoriented graphs, use BFS shortest cycle through edge (s,t)
    // to find a true face. The 3-element prefix {s,t,r} is unreliable
    // because next(t,s) depends on arbitrary neighbor order.
    outer_face = shortest_cycle({s,t},face_max);
  }

  return tutte_layout(outer_face);
}

vector<coord2d> PlanarGraphView::tutte_layout(const face_t& outer_face) const
{
  unsigned int Nface = outer_face.size();
  vector<coord2d> initial_coords(N);
  for(unsigned int i=0;i<Nface;i++){
    // the outer face is layed out CW:
    initial_coords[outer_face[i]] = coord2d(sin(i*2*M_PI/double(Nface)),cos(i*2*M_PI/double(Nface)));
  }

//  cout << "g = " << *this << endl;

// TODO: There seems to currently be a bug in the direct solver! When did this appear? Need to investigate and fix.
//  initial_coords = tutte_layout_direct(outer_face,initial_coords);
  return tutte_layout_iterative(outer_face,initial_coords);
}
#ifdef HAS_MKL
# include <mkl.h>
# include <mkl_pardiso.h>
# include <mkl_types.h>
#else
# include "../contrib/mgmres.hpp"
#endif

vector<coord2d> PlanarGraphView::tutte_layout_direct(const face_t& outer_face, const vector<coord2d>& initial_coords) const
{
  vector<coord2d> result(N);

  // Construct right hand side: rhs_i = {0,0} everywhere except for outer-face nodes with fixed solutions 
  double rhs[N*2], x[N*2];

  memset(rhs,0,N*2*sizeof(double));
  for(int i=0;i<N;i++){
    rhs[i]   = initial_coords[i].first;
    rhs[N+i] = initial_coords[i].second;
  }
  memcpy(x,rhs,N*2*sizeof(double));

  // Construct matrix I-1/3*A in CSR form, where A is adjacency, except for rows A_i=\delta_{ij} for i an outer-face node
  double A[4*N];
  int IA[N+1], JA[4*N];
  int nz = 0;
  {
    double *Afull = (double*)calloc(N*N,sizeof(double));
    assert(Afull != 0);
    memset(Afull,0,N*N*sizeof(double));
    for(node_t u=0;u<N;u++){
      Afull[u*(N+1)] = 1.0;
      for(int i=0;i<(*this)[u].size();i++){
        const node_t& v((*this)[u][i]);
        Afull[u*N+v] = -1.0L/3.0;
      }
    }
    for(int i=0;i<outer_face.size();i++){
      const node_t& u(outer_face[i]);
      for(node_t v=0;v<N;v++) Afull[u*N+v] = (u==v)? 1.0 : 0.0;
    }

    // Write in Compressed Sparse Row format
    for(node_t u=0;u<N;u++){
      IA[u] = nz;
      for(node_t v=0;v<N;v++)
        if(Afull[u*N+v]!=0){
          A[nz] = Afull[u*N+v];
          JA[nz] = v;
          nz++;
        }
    }
    free(Afull);
  }
  IA[N] = nz;

  // Solve sparse linear system for x-coordinates and y-coordinates 
#ifdef HAS_MKL
  {
    /* Matrix data. */
    MKL_INT n = N;
    MKL_INT mtype = 11;		/* Real unsymmetric matrix */
    MKL_INT nrhs = 2;		/* Number of right hand sides. */
    void *pt[64];        /* Internal solver memory pointer pt, */
    MKL_INT iparm[64];  /* Pardiso control parameters. */
    MKL_INT maxfct = 1, mnum = 1, phase, error = 0, msglvl = 0;
    /* Auxiliary variables. */
    double ddum;			/* Double dummy */
    MKL_INT idum;			/* Integer dummy. */
    /* -------------------------------------------------------------------- */
    /* .. Setup Pardiso control parameters. */
    /* -------------------------------------------------------------------- */
    memset(iparm,0,64*sizeof(MKL_INT));
    memset(pt,0,64*sizeof(void*));

    iparm[0] = 1;			/* No solver default */
    iparm[1] = 2;			/* Fill-in reordering from METIS */
    /* Numbers of processors, value of OMP_NUM_THREADS */
    iparm[2] = 0;
    iparm[3] = 1;
    iparm[7] = 200;		/* Max numbers of iterative refinement steps */
    iparm[9] = 13;		/* Perturb the pivot elements with 1E-13 */
    iparm[10] = 1;		/* Use nonsymmetric permutation and scaling MPS */
    iparm[12] = 1;		/* Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
    iparm[17] = -1;		/* Output: Number of nonzeros in the factor LU */
    iparm[18] = -1;		/* Output: Mflops for LU factorization */
    iparm[19] = 0;		/* Output: Numbers of CG Iterations */
    iparm[34] = 1;		/* PARDISO use C-style indexing for ia and ja arrays */

    /* -------------------------------------------------------------------- */
    /* .. Reordering and Symbolic Factorization. This step also allocates */
    /* all memory that is necessary for the factorization. */
    /* -------------------------------------------------------------------- */
    phase = 11;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
	     &n, A, IA, JA, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if (error != 0)
      {
	printf ("\nERROR during symbolic factorization of Tutte embedding: %d", error);
	exit (1);
      }
    printf ("\nReordering completed ... ");
    printf ("\nNumber of nonzeros in factors = %d", iparm[17]);
    printf ("\nNumber of factorization MFLOPS = %d", iparm[18]);
    /* -------------------------------------------------------------------- */
    /* .. Numerical factorization. */
    /* -------------------------------------------------------------------- */
    phase = 22;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
	     &n, A, IA, JA, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if (error != 0)
      {
	printf ("\nERROR during numerical factorization of Tutte embedding: %d", error);
	exit (2);
      }
    printf ("\nFactorization completed ... ");
    /* -------------------------------------------------------------------- */
    /* .. Back substitution and iterative refinement. */
    /* -------------------------------------------------------------------- */
    phase = 33;

    //    printf ("\n\nSolving system...\n");
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
	     &n, A, IA, JA, &idum, &nrhs, iparm, &msglvl, rhs, x, &error);
    if (error != 0)
      {
	printf ("\nERROR during direct sparse solution of Tutte embedding: %d", error);
	exit (3);
      }


    /* -------------------------------------------------------------------- */
    /* .. Termination and release of memory. */
    /* -------------------------------------------------------------------- */
    phase = -1;			/* Release internal memory. */
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
	     &n, &ddum, IA, JA, &idum, &nrhs,
	     iparm, &msglvl, &ddum, &ddum, &error);
  }

#else
  pmgmres_ilu_cr(N,nz,IA,JA,A,x,  rhs,  50000,N-1,1e-13,1e-13);
  pmgmres_ilu_cr(N,nz,IA,JA,A,x+N,rhs+N,50000,N-1,1e-13,1e-13);
#endif

  for(node_t u=0;u<N;u++) result[u] = coord2d(x[u],x[u+N]);

  return result;
}


vector<coord2d> PlanarGraphView::tutte_layout_iterative(const face_t& outer_face, const vector<coord2d>& initial_coords) const
{
  vector<coord2d> xys(initial_coords.begin(), initial_coords.end()), newxys(N);
  vector<bool> fixed(N);

  //  cerr << "tutte_layout: Outer face: " << outer_face << endl;

  unsigned int Nface = outer_face.size();
  for(unsigned int i=0;i<Nface;i++)
    fixed[outer_face[i]] = true;

  
  bool converged = false;
  const unsigned int TUTTE_MAX_ITERATION = 10000000;
  const double TUTTE_CONVERGENCE = 5e-4;
  unsigned int i;
  double max_change;
  for(i=0;!converged && i<TUTTE_MAX_ITERATION; i++){

    max_change = 0;

    for(node_t u=0;u<N;u++)
      if(fixed[u]){
        newxys[u] = xys[u];
      } else {
        auto ns = (*this)[u];
        coord2d neighbour_sum(0.0);
        
        for(int i=0;i<ns.size();i++) neighbour_sum += xys[ns[i]];
        newxys[u] = xys[u]*0.15 + (neighbour_sum/ns.size())*0.85;
        
        // Did the solution converge yet?
        double neighbour_dist = 0;
        for(size_t i=0;i<ns.size();i++) neighbour_dist += (xys[u]-xys[ns[i]]).norm()/ns.size();
        if(neighbour_dist > 0.0){ // let's not divide by zero
          double relative_change = (xys[u]-newxys[u]).norm()/neighbour_dist;
          if(relative_change > max_change) max_change = relative_change;
        }
      }
      
    if(max_change <= TUTTE_CONVERGENCE) converged = true;
    xys = newxys;
  }
  if(i>=TUTTE_MAX_ITERATION){
    printf("Planar Tutte embedding failed to converge. Increase TUTTE_MAX_ITERATION. ");
    cout << "layout = " << xys << ";\n";
    abort();
  }
  //  cerr << "Tutte layout of "<<N<<" vertices converged after " << i << " iterations, with maximal relative change " << max_change << endl;
  // Test that points are distinct
  ToleranceLess lt(0.0);
  set<coord2d,ToleranceLess> point_set(xys.begin(),xys.end(),lt);
  if(point_set.size() != N){
    fprintf(stderr,"Tutte layout failed: only %d unique coordinates out of %d vertices (up to tolerance %g).\n",
	    int(point_set.size()),N,0.0);
    abort();
  }
  return xys;
}

// spherical_projection moved to layout2d.cc


coord2d GraphView::centre2d(const vector<coord2d>& layout) const {
  coord2d centre(0,0);
  for(node_t u=0;u<layout.size();u++) centre += layout[u];
  return centre/double(layout.size());
}

coord3d GraphView::centre3d(std::span<const coord3d> layout) const {
  coord3d centre(0,0,0);
  for(node_t u=0;u<layout.size();u++) centre += layout[u];
  return centre/double(layout.size());
}


// to_latex and to_povray moved to layout2d.cc

