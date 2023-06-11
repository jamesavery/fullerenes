//Pentagons = 0
//Hexagons = 1
constexpr __constant__ float optimal_corner_cos_angles[2] = {-0.30901699437494734, -0.5}; 
constexpr __constant__ float optimal_bond_lengths[3] = {1.479, 1.458, 1.401}; 
constexpr __constant__ float optimal_dih_cos_angles[8] = {0.7946545571495363, 0.872903607049519, 0.872903607049519, 0.9410338472965512, 0.8162879359966257, 0.9139497166300941, 0.9139497166300941, 1.}; 

#if SEMINARIO_FORCE_CONSTANTS==1
constexpr __constant__ float angle_forces[2] = {207.924,216.787}; 
constexpr __constant__ float bond_forces[3] = {260.0, 353.377, 518.992}; 
constexpr __constant__ float dih_forces[4] = {35.0,65.0,3.772,270.0}; 
constexpr __constant__ float flat_forces[3] = {0., 0., 0.};
#else
constexpr __constant__ float angle_forces[2] = {100.0,100.0}; 
constexpr __constant__ float bond_forces[3] = {260.0,390.0,450.0}; 
constexpr __constant__ float dih_forces[4] = {35.0,65.0,85.0,270.0}; 
constexpr __constant__ float flat_forces[3] = {0., 0., 0.};
#endif

struct Constants{
    #if USE_CONSTANT_INDICES
    uchar4 i_f_bond;
    uchar4 i_f_inner_angle;
    uchar4 i_f_inner_dihedral;
    uchar4 i_f_outer_angle_m;
    uchar4 i_f_outer_angle_p;
    uchar4 i_f_outer_dihedral;
    uchar4 i_r0;
    uchar4 i_angle0;
    uchar4 i_outer_angle_m0;
    uchar4 i_outer_angle_p0;
    uchar4 i_inner_dih0;
    uchar4 i_outer_dih0_a;
    uchar4 i_outer_dih0_m;
    uchar4 i_outer_dih0_p;

    //Load force constants from neighbouring face information.
    constexpr INLINE device_real_t r0(const uint8_t j) const {return  optimal_bond_lengths[d_get(i_f_bond, j)];}
    constexpr INLINE device_real_t angle0(const uint8_t j)  const {return  optimal_corner_cos_angles[d_get(i_f_inner_angle, j)];}
    constexpr INLINE device_real_t inner_dih0(const uint8_t j) const {return  optimal_dih_cos_angles[d_get(i_f_inner_dihedral, j)];}
    constexpr INLINE device_real_t outer_angle_m0(const uint8_t j) const {return  optimal_corner_cos_angles[d_get(i_f_outer_angle_m, j)];}
    constexpr INLINE device_real_t outer_angle_p0(const uint8_t j) const {return  optimal_corner_cos_angles[d_get(i_f_outer_angle_p, j)];}
    constexpr INLINE device_real_t outer_dih0_a(const uint8_t j) const {return  optimal_dih_cos_angles[d_get(i_f_outer_dihedral, j)];}
    constexpr INLINE device_real_t outer_dih0_m(const uint8_t j) const {return  optimal_dih_cos_angles[d_get(i_r0, j)];}
    constexpr INLINE device_real_t outer_dih0_p(const uint8_t j) const {return  optimal_dih_cos_angles[d_get(i_angle0, j)];}
    constexpr INLINE device_real_t f_bond(const uint8_t j) const {return  bond_forces[d_get(i_outer_angle_m0, j)];}
    constexpr INLINE device_real_t f_inner_angle(const uint8_t j) const {return  angle_forces[d_get(i_outer_angle_p0, j)];}
    constexpr INLINE device_real_t f_inner_dihedral(const uint8_t j) const {return  dih_forces[d_get(i_inner_dih0, j)];}
    constexpr INLINE device_real_t f_outer_angle_m(const uint8_t j) const {return  angle_forces[d_get(i_outer_dih0_a, j)];}
    constexpr INLINE device_real_t f_outer_angle_p(const uint8_t j) const {return  angle_forces[d_get(i_outer_dih0_m, j)];}
    constexpr INLINE device_real_t f_outer_dihedral(const uint8_t j) const {return  dih_forces[d_get(i_outer_dih0_p, j)];}
    constexpr INLINE device_real_t f_flat() const {return 5e2;}
    #else
    device_coord3d f_bond;
    device_coord3d f_inner_angle;
    device_coord3d f_inner_dihedral;
    device_coord3d f_outer_angle_m;
    device_coord3d f_outer_angle_p;
    device_coord3d f_outer_dihedral;
    device_real_t f_flat = 1e2;
    
    device_coord3d r0;
    device_coord3d angle0;
    device_coord3d outer_angle_m0;
    device_coord3d outer_angle_p0;
    device_coord3d inner_dih0;
    device_coord3d outer_dih0_a;
    device_coord3d outer_dih0_m;
    device_coord3d outer_dih0_p;
    #endif

    __device__ __host__ __forceinline__ uint8_t face_index(uint8_t f1, uint8_t f2, uint8_t f3){
        return f1*4 + f2*2 + f3;
    }

    /**
     * @brief Constructor for the Constants struct
     *
     * @param G The IsomerBatch in which the graph information is read from
     * @param isomer_idx The index of the isomer that the current thread is a part of
     * @return Forcefield constants for the current node in the isomer_idx^th isomer in G
     */
    INLINE Constants(const IsomerBatch& G, const uint32_t isomer_idx){
        //Set pointers to start of fullerene.
        const DeviceCubicGraph FG(&G.cubic_neighbours[isomer_idx*blockDim.x*3]);
        device_node3 cubic_neighbours = {FG.cubic_neighbours[threadIdx.x*3], FG.cubic_neighbours[threadIdx.x*3 + 1], FG.cubic_neighbours[threadIdx.x*3 + 2]};
        //       m    p
        //    f5_|   |_f4
        //   p   c    b  m
        //       \f1/
        //     f2 a f3
        //        |
        //        d
        //      m/\p
        //       f6
        
        for (uint8_t j = 0; j < 3; j++) {
            //Faces to the right of arcs ab, ac and ad.
            
            uint8_t F1 = FG.face_size(threadIdx.x, d_get(cubic_neighbours, j)) - 5;
            uint8_t F2 = FG.face_size(threadIdx.x, d_get(cubic_neighbours, (j+1)%3)) -5;
            uint8_t F3 = FG.face_size(threadIdx.x, d_get(cubic_neighbours, (j+2)%3)) -5;
            
            //The faces to the right of the arcs ab, bm and bp in no particular order, from this we can deduce F4.
            uint8_t neighbour_F1 = FG.face_size(d_get(cubic_neighbours, j), FG.cubic_neighbours[d_get(cubic_neighbours, j)*3] ) -5;
            uint8_t neighbour_F2 = FG.face_size(d_get(cubic_neighbours, j), FG.cubic_neighbours[d_get(cubic_neighbours, j)*3 + 1] ) -5;
            uint8_t neighbour_F3 = FG.face_size(d_get(cubic_neighbours, j), FG.cubic_neighbours[d_get(cubic_neighbours, j)*3 + 2] ) -5;

            uint8_t F4 = neighbour_F1 + neighbour_F2 + neighbour_F3 - F1 - F3 ;
            
            //Load equillibirium distance, angles and dihedral angles from face information.
            #if USE_CONSTANT_INDICES
            d_set(i_r0,               j,  F3 + F1);
            d_set(i_angle0,           j,  F1);
            d_set(i_inner_dih0,       j,  face_index(F1, F2 , F3));
            d_set(i_outer_angle_m0,   j,  F3);
            d_set(i_outer_angle_p0,   j,  F1);
            d_set(i_outer_dih0_a,     j,  face_index(F3, F4, F1));
            d_set(i_outer_dih0_m,     j,  face_index(F4, F1, F3));
            d_set(i_outer_dih0_p,     j,  face_index(F1, F3, F4));
            
            //Load force constants from neighbouring face information.
            d_set(i_f_bond,           j,  bond_forces[F3 + F1]);
            d_set(i_f_inner_angle,    j,  angle_forces[F1]);
            d_set(i_f_inner_dihedral, j,  dih_forces[F1 + F2 + F3]);
            d_set(i_f_outer_angle_m,  j,  angle_forces[F3]);
            d_set(i_f_outer_angle_p,  j,  angle_forces[F1]);
            d_set(i_f_outer_dihedral, j,  dih_forces[F1 + F3 + F4]);
            #else
            d_set(r0,               j,  optimal_bond_lengths[F3 + F1]);
            d_set(angle0,           j,  optimal_corner_cos_angles[F1]);
            d_set(inner_dih0,       j,  optimal_dih_cos_angles[face_index(F1, F2 , F3)]);
            d_set(outer_angle_m0,   j,  optimal_corner_cos_angles[F3]);
            d_set(outer_angle_p0,   j,  optimal_corner_cos_angles[F1]);
            d_set(outer_dih0_a,     j,  optimal_dih_cos_angles[face_index(F3, F4, F1)]);
            d_set(outer_dih0_m,     j,  optimal_dih_cos_angles[face_index(F4, F1, F3)]);
            d_set(outer_dih0_p,     j,  optimal_dih_cos_angles[face_index(F1, F3, F4)]);
            
            //Load force constants from neighbouring face information.
            d_set(f_bond,           j,  bond_forces[F3 + F1]);
            d_set(f_inner_angle,    j,  angle_forces[F1]);
            d_set(f_inner_dihedral, j,  dih_forces[F1 + F2 + F3]);
            d_set(f_outer_angle_m,  j,  angle_forces[F3]);
            d_set(f_outer_angle_p,  j,  angle_forces[F1]);
            d_set(f_outer_dihedral, j,  dih_forces[F1 + F3 + F4]);
            #endif
        }
    }   
};