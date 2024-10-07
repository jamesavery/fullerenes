//Pentagons = 0
//Hexagons = 1


/* #if SEMINARIO_FORCE_CONSTANTS==1
constexpr float angle_forces[2] = {207.924,216.787}; 
constexpr float bond_forces[3] = {260.0, 353.377, 518.992}; 
constexpr float dih_forces[4] = {35.0,65.0,3.772,270.0}; 
constexpr float flat_forces[3] = {0., 0., 0.};
#else
#endif */


template <typename T, typename K>
struct Constants{
    
    TEMPLATE_TYPEDEFS(T,K);
    
    coord3d f_bond;
    coord3d f_inner_angle;
    coord3d f_inner_dihedral;
    coord3d f_outer_angle_m;
    coord3d f_outer_angle_p;
    coord3d f_outer_dihedral;
    real_t f_flat = 2e2;
    
    coord3d r0;
    coord3d angle0;
    coord3d outer_angle_m0;
    coord3d outer_angle_p0;
    coord3d inner_dih0;
    coord3d outer_dih0_a;
    coord3d outer_dih0_m;
    coord3d outer_dih0_p;

    /**
     * @brief Constructor for the Constants struct
     *
     * @param G The IsomerBatch in which the graph information is read from
     * @param isomer_idx The index of the isomer that the current thread is a part of
     * @return Forcefield constants for the current node in the isomer_idx^th isomer in G
     */

    inline Constants(const Span<std::array<K,3>> cubic_neighbours, K idx){

        constexpr real_t optimal_corner_cos_angles[2] = {-0.30901699437494734, -0.5}; 
        constexpr real_t optimal_bond_lengths[3] = {1.479, 1.458, 1.401}; 
        constexpr real_t optimal_dih_cos_angles[8] = {0.7946545571495363, 0.872903607049519, 0.872903607049519, 0.9410338472965512, 0.8162879359966257, 0.9139497166300941, 0.9139497166300941, 1.}; 
        constexpr real_t angle_forces[2] = {100.0,100.0}; 
        constexpr real_t bond_forces[3] = {260.0,390.0,450.0}; 
        constexpr real_t dih_forces[4] = {35.0,65.0,85.0,270.0}; 
        //constexpr float flat_forces[3] = {0., 0., 0.};
        //Set pointers to start of fullerene.

        auto face_index = [&](int f1, int f2, int f3){
            return f1*4 + f2*2 + f3;
        };


        const DeviceCubicGraph<K> FG(cubic_neighbours);
        node3 neighbours = {FG[idx][0], FG[idx][1], FG[idx][2]};
        //       m    p
        //    f5_|   |_f4
        //   p   c    b  m
        //       \f1/
        //     f2 a f3
        //        |
        //        d
        //      m/\p
        //       f6
        
        for (int j = 0; j < 3; j++) {
            //Faces to the right of arcs ab, ac and ad.
            
            int F1 = FG.face_size(idx, neighbours[j]) - 5;
            int F2 = FG.face_size(idx, neighbours[(j+1)%3]) -5;
            int F3 = FG.face_size(idx, neighbours[(j+2)%3]) -5;
            
            //The faces to the right of the arcs ab, bm and bp in no particular order, from this we can deduce F4.
            int neighbour_F1 = FG.face_size(neighbours[j], FG[neighbours[j]][0] ) -5;
            int neighbour_F2 = FG.face_size(neighbours[j], FG[neighbours[j]][1] ) -5;
            int neighbour_F3 = FG.face_size(neighbours[j], FG[neighbours[j]][2] ) -5;

            int F4 = (neighbour_F1 + neighbour_F2 + neighbour_F3 - F1 - F3) ;
            
            //Load equillibirium distance, angles and dihedral angles from face information.

            r0[j]             =  (T)optimal_bond_lengths[F3 + F1];
            angle0[j]         =  (T)optimal_corner_cos_angles[F1];
            inner_dih0[j]     =  (T)optimal_dih_cos_angles[face_index(F1, F2 , F3)];
            outer_angle_m0[j] =  (T)optimal_corner_cos_angles[F3];
            outer_angle_p0[j] =  (T)optimal_corner_cos_angles[F1];

            //outer_dih0_a[j]   =  (T)optimal_dih_cos_angles[(neighbour_F1 -F1 - F3)];
            outer_dih0_a[j]   =  (T)optimal_dih_cos_angles[face_index(F3, F4, F1)];
            outer_dih0_m[j]   =  (T)optimal_dih_cos_angles[face_index(F4, F1, F3)];
            outer_dih0_p[j]   =  (T)optimal_dih_cos_angles[face_index(F1, F3, F4)];
            
            //Load force constants from neighbouring face information.
            f_bond[j]           =  (T)bond_forces[F3 + F1];
            f_inner_angle[j]    =  (T)angle_forces[F1];
            f_inner_dihedral[j] =  (T)dih_forces[F1 + F2 + F3];
            f_outer_angle_m[j]  =  (T)angle_forces[F3];
            f_outer_angle_p[j]  =  (T)angle_forces[F1];
            f_outer_dihedral[j] =  (T)dih_forces[F1 + F3 + F4];
            /* 
            */
        } 
    }   
};
