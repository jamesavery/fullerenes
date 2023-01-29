#include <cuda.h>
#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>
#include <cooperative_groups/scan.h>
#include "cuda_runtime.h"
#include <cuda_runtime_api.h>
#include "fullerenes/gpu/cuda_definitions.h"
#include "fullerenes/gpu/cu_array.hh"
#include "fullerenes/gpu/isomer_batch.hh"
#include "fullerenes/gpu/kernels.hh"

namespace gpu_kernels{
namespace isomerspace_forcefield{
#include "device_includes.cu"

// This struct was made to reduce signature cluttering of device functions, it is simply a container for default arguments which are shared between functions
template <ForcefieldType T>
struct ForceField{
    DEVICE_TYPEDEFS
    
    const NodeNeighbours node_graph;         //Contains face-information and neighbour-information. Both of which are constant in the lifespan of this struct. 
    const Constants constants;          //Contains force-constants and equillibrium-parameters. Constant in the lifespan of this struct.

    size_t node_id = threadIdx.x;
    real_t* sdata;                      //Pointer to start of L1 cache array, used exclusively for reduction.

    __device__ ForceField(  const NodeNeighbours &G,
                            const Constants &c, 
                            real_t* sdata): node_graph(G), constants(c), sdata(sdata) {}



struct FaceData{
    coord3d Xa;
    symMat3 A;
    coord3d n_f; //normalized normal vector to face-plane
    real_t lambda_f; //Smallest eigenvalue defining the flatness of the face
    coord3d lambdas;
    coord3d centroid;
    device_node3 face_neighbours;
    //84 + 107 FLOPS
    INLINE FaceData(const coord3d* X, const NodeNeighbours& G){
        face_neighbours = G.face_neighbours;
        Xa = X[threadIdx.x];
        //There are only blockDim.x/2 + 2 faces. (Nf  =  N/2 + 1)
        if(threadIdx.x < blockDim.x/2 + 2){
            coord3d Xf[6] = {X[G.face_nodes[0]], X[G.face_nodes[1]] , X[G.face_nodes[2]] , X[G.face_nodes[3]] , X[G.face_nodes[4]] };
            //If pentagon set to 0 otherwise get the 6th node coordinates.
            if(G.face_size == 6){Xf[5] = X[G.face_nodes[5]];} else {Xf[5] = {(real_t)0., (real_t)0., (real_t)0.};}
            centroid = (Xf[0] + Xf[1] + Xf[2] + Xf[3] + Xf[4] + Xf[5]) / (device_real_t)G.face_size;
            //Centralise coordinate system to centroid of the face
            Xf[0] -= centroid; Xf[1] -= centroid; Xf[2] -= centroid; Xf[3] -= centroid; Xf[4] -= centroid;  
            if(G.face_size == 6){Xf[5] -= centroid;}
            auto a = Xf[0].x * Xf[0].x + Xf[1].x * Xf[1].x + Xf[2].x * Xf[2].x + Xf[3].x * Xf[3].x + Xf[4].x * Xf[4].x + Xf[5].x * Xf[5].x,
                 b = Xf[0].x * Xf[0].y + Xf[1].x * Xf[1].y + Xf[2].x * Xf[2].y + Xf[3].x * Xf[3].y + Xf[4].x * Xf[4].y + Xf[5].x * Xf[5].y,
                 c = Xf[0].x * Xf[0].z + Xf[1].x * Xf[1].z + Xf[2].x * Xf[2].z + Xf[3].x * Xf[3].z + Xf[4].x * Xf[4].z + Xf[5].x * Xf[5].z,
                 d = Xf[0].y * Xf[0].y + Xf[1].y * Xf[1].y + Xf[2].y * Xf[2].y + Xf[3].y * Xf[3].y + Xf[4].y * Xf[4].y + Xf[5].y * Xf[5].y,
                 e = Xf[0].y * Xf[0].z + Xf[1].y * Xf[1].z + Xf[2].y * Xf[2].z + Xf[3].y * Xf[3].z + Xf[4].y * Xf[4].z + Xf[5].y * Xf[5].z,
                 f = Xf[0].z * Xf[0].z + Xf[1].z * Xf[1].z + Xf[2].z * Xf[2].z + Xf[3].z * Xf[3].z + Xf[4].z * Xf[4].z + Xf[5].z * Xf[5].z;
            //Xf * Xf^T In closed form.
            A = symMat3(a,b,c,d,e,f);

            //A is positive-semi-definite so all eigenvalues are non-negative
            lambdas = A.eigenvalues();
            lambda_f = d_min(d_min(lambdas.x, lambdas.y), lambdas.z);  
        }
    }
    //3 FLOPs
        /**
     * Computes the harmonic energy contribution of one term.
     *
     * @param[in] p0 Equillibrium parameter
     * @param[in] p Current parameter
     * @return Hooke's law harmonic energy contribution of the term
     */
    INLINE real_t harmonic_energy(const real_t p0, const real_t p) const{
        return (real_t)0.5*(p-p0)*(p-p0);
    }


    /** @brief Compute the flatness of the threadIdx^th face in the isomer
     *  @return The flatness of the threadIdx^th face in the isomer
     */
    INLINE real_t flatness() const {return threadIdx.x < blockDim.x/2 + 2 ? lambda_f : (real_t)0.;}

    //4 FLOPs
    /**
     * Computes the harmonic energy gradient contribution of one term.
     *
     * @param[in] p0 Equillibrium parameter
     * @param[in] p Current parameter
     * @param[in] gradp Gradient of the parameter w.r.t. the particle position
     * @return Hooke's law harmonic energy gradient contribution of the term
     */
    INLINE coord3d  harmonic_energy_gradient(const real_t p0, const real_t p, const coord3d gradp) const{
        return (p-p0)*gradp;     
    }


    /**
 * @brief Compute the flatness energy contribution of the threadIdx^th face in the isomer.
 *
 * @param c The forcefield constants for the threadIdx^th node.
 * @return The flatness energy.
 */
    INLINE real_t flatness_energy(const Constants& c) const {
        #if USE_CONSTANT_INDICES
        return c.f_flat() * harmonic_energy(flatness(),(real_t)0.);
        #else
        return c.f_flat * harmonic_energy(flatness(),(real_t)0.);
        #endif
    }

    /**
     * @brief Compute the gradient of the flatness w.r.t to the threadIdx^th atom in the isomer.
     * @param c The forcefield constants for the threadIdx^th node.
     * @param cache A pointer to a cache of minimum size Nf * 2 * sizeof(coord3d) bytes.
     * @return The flatness energy gradient.
     */
    INLINE coord3d flatness_gradient(const Constants& c, coord3d* cache) const {
        coord3d* centroids = reinterpret_cast<coord3d*>(cache);
        coord3d* norms = reinterpret_cast<coord3d*>(cache + blockDim.x/2 + 2);
        if(threadIdx.x < blockDim.x/2 + 2){
            centroids[threadIdx.x] = centroid;
            norms[threadIdx.x] = A.eigenvector(lambda_f);
        }
        BLOCK_SYNC

        coord3d grad = {(real_t)0., (real_t)0., (real_t)0.};
        for(unsigned char j = 0; j < 3; j++) grad += dot(Xa - centroids[d_get(face_neighbours,j)], norms[d_get(face_neighbours,j)]) * norms[d_get(face_neighbours,j)];
        #if USE_CONSTANT_INDICES
        return c.f_flat() * (real_t)2. * grad;
        #else
        return c.f_flat * (real_t)2. * grad;
        #endif

    }
};

//Container for all energy and gradient evaluations with respect to an arc, eg. AB, AC or AD.
struct ArcData{
    //124 FLOPs;
    uint8_t j;
    /**
     * @brief Construct a new ArcData object
     * @param j The index of the arc, eg. 0 for ab, 1 for ac and 2 for ad.
     * @param X The coordinates of all nodes in the isomer.
     * @param G The neighbour information for the threadIdx^th node.
     * @return A new ArcData object.
    */
    INLINE ArcData(const uint8_t j, const coord3d* __restrict__ X, const NodeNeighbours& G){  
        __builtin_assume(j < 3); 
        this->j = j;   
        node_t a = threadIdx.x;
        real_t r_rmp;
        coord3d ap, am, ab, ac, ad, mp;
        coord3d X_a = X[a], X_b = X[d_get(G.cubic_neighbours,j)];

        //Compute the arcs ab, ac, ad, bp, bm, ap, am, mp, bc and cd
        ab = (X_b - X_a);  r_rab = bond_length(ab); ab_hat = r_rab * ab;
        ac = (X[d_get(G.cubic_neighbours,(j+1)%3)] - X_a); r_rac = bond_length(ac); ac_hat = r_rac * ac; rab = non_resciprocal_bond_length(ab);
        ad = (X[d_get(G.cubic_neighbours,(j+2)%3)] - X_a); r_rad = bond_length(ad); ad_hat = r_rad * ad;
        
        coord3d bp = (X[d_get(G.next_on_face,j)] - X_b); bp_hat = unit_vector(bp);
        coord3d bm = (X[d_get(G.prev_on_face,j)] - X_b); bm_hat = unit_vector(bm);

        ap = bp + ab; r_rap = bond_length(ap); ap_hat = r_rap * ap;
        am = bm + ab; r_ram = bond_length(am); am_hat = r_ram * am;
        mp = bp - bm; r_rmp = bond_length(mp); mp_hat = r_rmp * mp;

        bc_hat = unit_vector(ac - ab);
        cd_hat = unit_vector(ad - ac);

        //Compute inverses of some arcs, these are subject to be omitted if the equations are adapted appropriately with inversion of signs.
        ba_hat = -ab_hat;
        mb_hat = -bm_hat;
        pa_hat = -ap_hat;
        pb_hat = -bp_hat;
    }

    //3 FLOPs
    /**
     * @brief Compute the harmonic energy contribution from one parameter.
     * @param p0 The equillibrium value of the parameter.
     * @param p The current value of the parameter.
    */
    INLINE real_t harmonic_energy(const real_t p0, const real_t p) const{
        return (real_t)0.5*(p-p0)*(p-p0);
    }
    //4 FLOPs
    /**
     * @brief Compute the harmonic energy gradient contribution from one parameter.
     * @param p0 The equillibrium value of the parameter.
     * @param p The current value of the parameter.
     * @param gradp The gradient of the parameter with respect to the node position.
    */
    INLINE coord3d  harmonic_energy_gradient(const real_t p0, const real_t p, const coord3d gradp) const{
        return (p-p0)*gradp;     
    }
    //1 FLOP
    /**
     * @brief Compute the bond length of the main arc ab or ac or ad. for j = 0, 1 or 2 respectively.
     * @return The bond length.
    */
    INLINE real_t bond() const {return rab;}
    //5 FLOPs
    /**
     * @brief Compute the cosine of the angle between the main arc and the next arc, (ab,ac), (ac,ad), (ad,ab). For j = 0, 1 or 2 respectively.
     * @return The cosine of the angle.
    */
    INLINE real_t angle() const {return dot(ab_hat,ac_hat);}

    INLINE real_t normalized_angle_err() const {return acos(dot(ab_hat,ac_hat));}

    //Returns outer angle m, used only diagnostically.
    INLINE real_t outer_angle_m() const {return -dot(ab_hat, bm_hat);} //Compute outer angle. ab,bm

    //Returns outer angle p, used only diagnostically.
    INLINE real_t outer_angle_p() const{return -dot(ab_hat, bp_hat);} //Compute outer angle. ab,bp

    //Returns the inner dihedral angle for the current arc. Used here only for energy calculation, 
    //otherwise embedded in dihedral computation because the planes and angles that make up the dihedral angle computation are required for derivative computation.
    //50 FLOPs
    /**
     * @brief Compute the dihedral angle between the faces (abc,bcd), (acd,bcd) and (abd,bcd). For j = 0, 1 or 2 respectively. 
     * @return The dihedral angle.
    */
    INLINE real_t dihedral() const 
    { 
        coord3d nabc, nbcd; real_t cos_b, cos_c, r_sin_b, r_sin_c;
        cos_b = dot(ba_hat,bc_hat);  r_sin_b = (real_t)1.0/sqrt((real_t)1.0 - cos_b*cos_b); nabc = cross(ba_hat, bc_hat) * r_sin_b;
        cos_c = dot(-bc_hat,cd_hat); r_sin_c = (real_t)1.0/sqrt((real_t)1.0 - cos_c*cos_c); nbcd = cross(-bc_hat,cd_hat) * r_sin_c;
        return dot(nabc, nbcd);
    }
    //Returns the Outer-dihedral-a wrt. current arc, only accessed diagnostically (internal coordinate).
    /**
     * @brief Compute the dihedral angle between the faces $(b-a-b_m, a-b_m-b_p)$, $(c-a-c_m, a-c_m-c_p)$ and $(d-a-d_m, a-d_m-d_p)$. For j = 0, 1 or 2 respectively.
     * @return The dihedral angle.
    */
    INLINE real_t outer_dihedral_a() const
    {
        coord3d nbam_hat, namp_hat; real_t cos_a, cos_m, r_sin_a, r_sin_m;
        cos_a = dot(ab_hat,am_hat); r_sin_a = (real_t)1.0/sqrt((real_t)1.0 - cos_a*cos_a); nbam_hat = cross(ab_hat,am_hat) * r_sin_a;
        cos_m = dot(-am_hat,mp_hat); r_sin_m = (real_t)1.0/sqrt((real_t)1.0 - cos_m*cos_m); namp_hat = cross(-am_hat,mp_hat) * r_sin_m;
        real_t cos_beta = dot(nbam_hat, namp_hat); //Outer Dihedral angle bam, amp
        return cos_beta;
    }
    //Returns the Outer-dihedral-m wrt. current arc, only accessed diagnostically (internal coordinate).
    /**
     * @brief Compute the dihedral angle between the faces $(b-b_m-b_p, b_m-b_p-a)$, $(c-c_m-c_p, c_m-c_p-a)$ and $(d-d_m-d_p, d_m-d_p-a)$. For j = 0, 1 or 2 respectively.
     * @return The dihedral angle.
    */
    INLINE real_t outer_dihedral_m() const
    {
        coord3d nbmp_hat, nmpa_hat; real_t cos_m, cos_p, r_sin_m, r_sin_p;
        cos_m = dot(mb_hat,mp_hat);  r_sin_m = (real_t)1.0/sqrt((real_t)1.0 - cos_m*cos_m); nbmp_hat = cross(mb_hat,mp_hat) * r_sin_m;
        cos_p = dot(-mp_hat,pa_hat); r_sin_p = (real_t)1.0/sqrt((real_t)1.0 - cos_p*cos_p); nmpa_hat = cross(-mp_hat,pa_hat) * r_sin_p;
        //Cosine to the outer dihedral angle constituted by the planes bmp and mpa
        real_t cos_beta = dot(nbmp_hat, nmpa_hat); //Outer dihedral angle bmp,mpa.
        return cos_beta;    
    }
    //Returns the Outer-dihedral-p wrt. current arc, only accessed diagnostically (internal coordinate).
    /**
     * @brief Compute the dihedral angle between the faces $(b-b_p-a, b_p-a-b_m)$, $(c-c_p-a, c_p-a-c_m)$ and $(d-d_p-a, d_p-a-d_m)$. For j = 0, 1 or 2 respectively.
     * @return The dihedral angle.
    */
    INLINE real_t outer_dihedral_p() const
    {
        coord3d nbpa_hat, npam_hat; real_t cos_p, cos_a, r_sin_p, r_sin_a;
        cos_a = dot(ap_hat,am_hat);  r_sin_a = (real_t)1.0/sqrt((real_t)1.0 - cos_a*cos_a); npam_hat = cross(ap_hat,am_hat)  * r_sin_a;
        cos_p = dot(pb_hat,-ap_hat); r_sin_p = (real_t)1.0/sqrt((real_t)1.0 - cos_p*cos_p); nbpa_hat = cross(pb_hat,-ap_hat) * r_sin_p;
        real_t cos_beta = dot(nbpa_hat, npam_hat); //Outer dihedral angle bpa, pam.
        //Eq. 33 multiplied by harmonic term.
        return cos_beta;
    }
    
    // Chain rule terms for angle calculation
    //Computes gradient related to bending term. ~24 FLOPs
    /**
     * @brief Compute the gradient of the bending term.
     * @param c The constants for the threadIdx^th node.
     * @return The gradient of the bending term.
    */
    INLINE coord3d inner_angle_gradient(const Constants& c) const
    {   
        real_t cos_angle = angle(); //Inner angle of arcs ab,ac.
        coord3d grad = cos_angle * (ab_hat * r_rab + ac_hat * r_rac) - ab_hat * r_rac - ac_hat* r_rab; //Derivative of inner angle: Eq. 21. 
        #if USE_CONSTANT_INDICES
        return c.f_inner_angle(j) * harmonic_energy_gradient(c.angle0(j), cos_angle, grad); //Harmonic Energy Gradient: Eq. 21. multiplied by harmonic term.
        #else
        return d_get(c.f_inner_angle,j) * harmonic_energy_gradient(d_get(c.angle0,j), cos_angle, grad); //Harmonic Energy Gradient: Eq. 21. multiplied by harmonic term.
        #endif
    }
    //Computes gradient related to bending of outer angles. ~20 FLOPs
    /**
     * @brief Compute the gradient of the outer angle-m term.
     * @param c The constants for the threadIdx^th node.
     * @return The gradient of the outer angle-m term.
    */
    INLINE coord3d outer_angle_gradient_m(const Constants& c) const
    {
        real_t cos_angle = -dot(ab_hat, bm_hat); //Compute outer angle. ab,bm
        coord3d grad = (bm_hat + ab_hat * cos_angle) * r_rab; //Derivative of outer angles Eq. 30. Buster Thesis
        #if USE_CONSTANT_INDICES
        return c.f_outer_angle_m(j) * harmonic_energy_gradient(c.outer_angle_m0(j),cos_angle,grad); //Harmonic Energy Gradient: Eq. 30 multiplied by harmonic term.
        #else
        return d_get(c.f_outer_angle_m,j) * harmonic_energy_gradient(d_get(c.outer_angle_m0,j),cos_angle,grad); //Harmonic Energy Gradient: Eq. 30 multiplied by harmonic term.
        #endif
    }

    /**
     * @brief Compute the gradient of the outer angle-p term.
     * @param c The constants for the threadIdx^th node.
     * @return The gradient of the outer angle-p term.
    */
    INLINE coord3d outer_angle_gradient_p(const Constants& c) const
    {   
        real_t cos_angle = -dot(ab_hat, bp_hat); //Compute outer angle. ab,bp
        coord3d grad = (bp_hat + ab_hat * cos_angle) * r_rab; //Derivative of outer angles Eq. 28. Buster Thesis
        #if USE_CONSTANT_INDICES
        return c.f_outer_angle_p(j) * harmonic_energy_gradient(c.outer_angle_p0(j),cos_angle,grad); //Harmonic Energy Gradient: Eq. 28 multiplied by harmonic term.
        #else
        return d_get(c.f_outer_angle_p,j) * harmonic_energy_gradient(d_get(c.outer_angle_p0,j),cos_angle,grad); //Harmonic Energy Gradient: Eq. 28 multiplied by harmonic term.
        #endif
    }
    // Chain rule terms for dihedral calculation
    //Computes gradient related to dihedral/out-of-plane term. ~75 FLOPs
    /**
     * @brief Compute the gradient of the inner dihedral term.
     * @param c The constants for the threadIdx^th node.
     * @return The gradient of the inner dihedral term.
    */
    INLINE coord3d inner_dihedral_gradient(const Constants& c) const
    {
        coord3d nabc, nbcd; real_t cos_b, cos_c, r_sin_b, r_sin_c;
        cos_b = dot(ba_hat,bc_hat); r_sin_b = (real_t)1.0/sqrt((real_t)1.0 - cos_b*cos_b); nabc = cross(ba_hat, bc_hat) * r_sin_b;
        cos_c = dot(-bc_hat,cd_hat); r_sin_c = (real_t)1.0/sqrt((real_t)1.0 - cos_c*cos_c); nbcd = cross(-bc_hat,cd_hat) * r_sin_c;

        real_t cos_beta = dot(nabc, nbcd); //Inner dihedral angle from planes abc,bcd.
        real_t cot_b = cos_b * r_sin_b * r_sin_b; //cos(b)/sin(b)^2

        //Derivative w.r.t. inner dihedral angle F and G in Eq. 26
        coord3d grad = cross(bc_hat, nbcd) * r_sin_b * r_rab - ba_hat * cos_beta * r_rab + (cot_b * cos_beta * r_rab) * (bc_hat - ba_hat * cos_b);
        #if USE_CONSTANT_INDICES
        return c.f_inner_dihedral(j) * harmonic_energy_gradient(c.inner_dih0(j), cos_beta, grad); //Eq. 26.
        #else
        return d_get(c.f_inner_dihedral,j) * harmonic_energy_gradient(d_get(c.inner_dih0,j), cos_beta, grad); //Eq. 26.
        #endif
    }

    //Computes gradient from dihedral angles constituted by the planes bam, amp ~162 FLOPs
    /**
     * @brief Compute the gradient of the outer dihedral-a term.
     * @param c The constants for the threadIdx^th node.
     * @return The gradient of the outer dihedral-a term.
    */
    INLINE coord3d outer_dihedral_gradient_a(const Constants& c) const
    {
        coord3d nbam_hat, namp_hat; real_t cos_a, cos_m, r_sin_a, r_sin_m;

        cos_a = dot(ab_hat,am_hat); r_sin_a = (real_t)1.0/sqrt((real_t)1.0 - cos_a*cos_a); nbam_hat = cross(ab_hat,am_hat) * r_sin_a;
        cos_m = dot(-am_hat,mp_hat); r_sin_m = (real_t)1.0/sqrt((real_t)1.0 - cos_m*cos_m); namp_hat = cross(-am_hat,mp_hat) * r_sin_m;
        
        real_t cos_beta = dot(nbam_hat, namp_hat); //Outer Dihedral angle bam, amp
        real_t cot_a = cos_a * r_sin_a * r_sin_a;
        real_t cot_m = cos_m * r_sin_m * r_sin_m;

        //Derivative w.r.t. outer dihedral angle, factorized version of Eq. 31.
        coord3d grad = cross(mp_hat,nbam_hat)*r_ram*r_sin_m - (cross(namp_hat,ab_hat)*r_ram + cross(am_hat,namp_hat)*r_rab)*r_sin_a +
                        cos_beta*(ab_hat*r_rab + r_ram * ((real_t)2.0*am_hat + cot_m*(mp_hat+cos_m*am_hat)) - cot_a*(r_ram*(ab_hat - am_hat*cos_a) + r_rab*(am_hat-ab_hat*cos_a)));
        
        //Eq. 31 multiplied by harmonic term.
        #if USE_CONSTANT_INDICES
        return c.f_outer_dihedral(j) * harmonic_energy_gradient(c.outer_dih0_a(j), cos_beta, grad);
        #else
        return d_get(c.f_outer_dihedral,j) * harmonic_energy_gradient(d_get(c.outer_dih0_a,j), cos_beta, grad);
        #endif
    }

    //Computes gradient from dihedral angles constituted by the planes nbmp, nmpa ~92 FLOPs
    /**
     * @brief Compute the gradient of the outer dihedral-m term.
     * @param c The constants for the threadIdx^th node.
     * @return The gradient of the outer dihedral-m term.
    */
    INLINE coord3d outer_dihedral_gradient_m(const Constants& c) const
    {
        coord3d nbmp_hat, nmpa_hat; real_t cos_m, cos_p, r_sin_m, r_sin_p;
        cos_m = dot(mb_hat,mp_hat);  r_sin_m = (real_t)1.0/sqrt((real_t)1.0 - cos_m*cos_m); nbmp_hat = cross(mb_hat,mp_hat) * r_sin_m;
        cos_p = dot(-mp_hat,pa_hat); r_sin_p = (real_t)1.0/sqrt((real_t)1.0 - cos_p*cos_p); nmpa_hat = cross(-mp_hat,pa_hat) * r_sin_p;
        
        //Cosine to the outer dihedral angle constituted by the planes bmp and mpa
        real_t cos_beta = dot(nbmp_hat, nmpa_hat); //Outer dihedral angle bmp,mpa.
        real_t cot_p = cos_p * r_sin_p * r_sin_p;
        
        //Derivative w.r.t. outer dihedral angle, factorized version of Eq. 32.
        coord3d grad = r_rap * (cot_p*cos_beta * (-mp_hat - pa_hat*cos_p) - cross(nbmp_hat, mp_hat)*r_sin_p - pa_hat*cos_beta );

        //Eq. 32 multiplied by harmonic term.
        #if USE_CONSTANT_INDICES
        return c.f_outer_dihedral(j) * harmonic_energy_gradient(c.outer_dih0_m(j), cos_beta, grad);
        #else
        return d_get(c.f_outer_dihedral,j) * harmonic_energy_gradient(d_get(c.outer_dih0_m,j), cos_beta, grad);
        #endif
    }

    //Computes gradient from dihedral angles constituted by the planes bpa, pam ~162 FLOPs
    /**
     * @brief Compute the gradient of the outer dihedral-p term.
     * @param c The constants for the threadIdx^th node.
     * @return The gradient of the outer dihedral-p term.
    */
    INLINE coord3d outer_dihedral_gradient_p(const Constants& c) const
    {
        coord3d nbpa_hat, npam_hat; real_t cos_p, cos_a, r_sin_p, r_sin_a;
        cos_a = dot(ap_hat,am_hat);  r_sin_a = (real_t)1.0/sqrt((real_t)1.0 - cos_a*cos_a); npam_hat = cross(ap_hat,am_hat)  * r_sin_a;
        cos_p = dot(pb_hat,-ap_hat); r_sin_p = (real_t)1.0/sqrt((real_t)1.0 - cos_p*cos_p); nbpa_hat = cross(pb_hat,-ap_hat) * r_sin_p;

        real_t cos_beta = dot(nbpa_hat, npam_hat); //Outer dihedral angle bpa, pam.
        real_t cot_p = cos_p * r_sin_p * r_sin_p;
        real_t cot_a = cos_a * r_sin_a * r_sin_a;

        //Derivative w.r.t. outer dihedral angle, factorized version of Eq. 33.
        coord3d grad = cross(npam_hat,pb_hat)*r_rap*r_sin_p - (cross(am_hat,nbpa_hat)*r_rap + cross(nbpa_hat,ap_hat)*r_ram)*r_sin_a +
                        cos_beta*(am_hat*r_ram + r_rap * ((real_t)2.0*ap_hat + cot_p*(pb_hat+cos_p*ap_hat)) - cot_a*(r_rap*(am_hat - ap_hat*cos_a) + r_ram*(ap_hat-am_hat*cos_a)));
        
        //Eq. 33 multiplied by harmonic term.
        #if USE_CONSTANT_INDICES
        return c.f_outer_dihedral(j) * harmonic_energy_gradient(c.outer_dih0_p(j), cos_beta, grad);
        #else
        return d_get(c.f_outer_dihedral,j) * harmonic_energy_gradient(d_get(c.outer_dih0_p,j), cos_beta, grad);
        #endif
    }

    // Internal coordinate gradients
    /**
     * @brief Compute the gradient of the bond length term.
     * @param c The constants for the threadIdx^th node.
     * @return The gradient of the bond length term.
    */
    INLINE coord3d bond_length_gradient(const Constants& c) const {
        #if USE_CONSTANT_INDICES
        return c.f_bond(j) * harmonic_energy_gradient(bond(),c.r0(j),ab_hat); 
        #else
        return d_get(c.f_bond,j) * harmonic_energy_gradient(bond(),d_get(c.r0,j),ab_hat); 
        #endif
    }
    //Sum of angular gradient components.
    /**
     * @brief Compute the sum of the gradients of the bending terms.
     * @param c The constants for the threadIdx^th node.
     * @return The sum of the gradients of the bending terms.
    */
    INLINE coord3d angle_gradient(const Constants& c) const { return inner_angle_gradient(c) + outer_angle_gradient_p(c) + outer_angle_gradient_m(c);}
    //Sum of inner and outer dihedral gradient components.
    /**
     * @brief Compute the sum of the gradients of the dihedral terms.
     * @param c The constants for the threadIdx^th node.
     * @return The sum of the gradients of the dihedral terms.
    */
    INLINE coord3d dihedral_gradient(const Constants& c) const { 
        switch (T)
        {
        case BUSTER:
            return inner_dihedral_gradient(c) + outer_dihedral_gradient_a(c) + outer_dihedral_gradient_m(c) + outer_dihedral_gradient_p(c);
        case WIRZ:
            return inner_dihedral_gradient(c);
        default:
            return inner_dihedral_gradient(c) + outer_dihedral_gradient_a(c) + outer_dihedral_gradient_m(c) + outer_dihedral_gradient_p(c);
        }
        return inner_dihedral_gradient(c) + outer_dihedral_gradient_a(c) + outer_dihedral_gradient_m(c) + outer_dihedral_gradient_p(c);
    }
    /**
     * @brief Compute the energy contribution of the bond length term.
     * @param c The constants for the threadIdx^th node.
     * @return The energy contribution of the bond length term.
    */
    INLINE real_t bond_energy(const Constants& c) const {
        #if USE_CONSTANT_INDICES
        return (real_t)0.5 *c.f_bond(j) *harmonic_energy(bond(),c.r0(j));
        #else
        return (real_t)0.5 *d_get(c.f_bond,j) *harmonic_energy(bond(),d_get(c.r0,j));
        #endif
    }
    /**
     * @brief Compute the total energy contribution of the bending terms.
     * @param c The constants for the threadIdx^th node.
     * @return The energy contribution of the bending terms.
    */
    INLINE real_t bend_energy(const Constants& c) const {
        #if USE_CONSTANT_INDICES
        return c.f_inner_angle(j)* harmonic_energy(angle(),c.angle0(j));
        #else
        return d_get(c.f_inner_angle,j)* harmonic_energy(angle(),d_get(c.angle0,j));
        #endif
    }

    /**
     * @brief Compute the total energy contribution of the dihedral terms.
     * @param c The constants for the threadIdx^th node.
     * @return The energy contribution of the dihedral terms.
    */
    INLINE real_t dihedral_energy(const Constants& c) const {
        #if USE_CONSTANT_INDICES
        return c.f_inner_dihedral(j)* harmonic_energy(dihedral(),c.inner_dih0(j));
        #else
        return d_get(c.f_inner_dihedral,j)* harmonic_energy(dihedral(),d_get(c.inner_dih0,j));
        #endif
    }
    //Harmonic energy contribution from bond stretching, angular bending and dihedral angle bending.
    //71 FLOPs
    /**
     * @brief Compute the total energy contribution of the bond length, bending and dihedral terms.
     * @param c The constants for the threadIdx^th node.
     * @return The energy contribution of the bond length, bending and dihedral terms.
    */
    INLINE real_t energy(const Constants& c) const {
        switch (T)
        {
        case FLAT_BOND:
            return bond_energy(c);
            break;
        default:
            return bond_energy(c) + bend_energy(c) + dihedral_energy(c);
            break;
        }
    }
    //Sum of bond, angular and dihedral gradient components.
    /**
     * @brief Compute the total gradient of the bond length, bending and dihedral terms w.r.t. the coordinates of the threadIdx^th node.
     * @param c The constants for the threadIdx^th node.
     * @return The gradient of the bond length, bending and dihedral terms w.r.t. the coordinates of the threadIdx^th node.
    */
    INLINE coord3d gradient(const Constants& c) const{
        switch (T)
        {
        case FLAT_BOND:
            return bond_length_gradient(c);
        default:
            return bond_length_gradient(c) + angle_gradient(c) + dihedral_gradient(c);
            break;
        }
        return bond_length_gradient(c) + angle_gradient(c) + dihedral_gradient(c);
    }


    //Reciprocal lengths of arcs ab, ac, am, ap.
    real_t
        rab,
        r_rab,
        r_rac,
        r_rad,
        r_ram,
        r_rap;

    //Base Arcs,
    coord3d
        ab,
        ac,
        ad;

    //All normalized arcs required to perform energy & gradient calculations.
    //Note that all these arcs are cyclical the arc ab becomes: ab->ac->ad,  the arc ac becomes: ac->ad->ab , the arc bc becomes: bc->cd->db (For iterations 0, 1, 2)
    //As such the naming convention here is related to the arcs as they are used in the 0th iteration.
    coord3d 
        ab_hat,
        ac_hat,
        ad_hat,
        bp_hat,
        bm_hat,
        am_hat,
        ap_hat,
        ba_hat,
        bc_hat,
        cd_hat,
        mp_hat,
        mb_hat,
        pa_hat,
        pb_hat;

    coord3d
        face_center, //Center of the face to the left of the arc a->b, a->b, a->c
        face_offset; //Difference between the node coordinates X and the face-center coordinates face_center
    
    coord3d A[3];
    
};


/**
 * @brief Compute the total gradient of the bond, flatness, bending and dihedral terms w.r.t. the coordinates of the threadIdx^th node.
 * @param c The constants for the threadIdx^th node.
 * @return The gradient of the bond, flatness, bending and dihedral terms w.r.t. the coordinates of the threadIdx^th node.
*/
INLINE coord3d gradient(coord3d* X) const {
    BLOCK_SYNC
    coord3d grad = {0.0, 0.0, 0.0};
    for (uint8_t j = 0; j < 3; j++ ){
        ArcData arc = ArcData(j, X, node_graph);
        grad += arc.gradient(constants);
    }
    switch (T)
    {
    case FLATNESS_ENABLED: {
        FaceData face(X, node_graph);
        auto face_grad = face.flatness_gradient(constants, reinterpret_cast<coord3d*>(sdata + Block_Size_Pow_2) + blockDim.x*2);
        return grad + face_grad;
        break;
        }
    case FLAT_BOND:{
        FaceData face(X, node_graph);
        auto face_grad = face.flatness_gradient(constants, reinterpret_cast<coord3d*>(sdata + Block_Size_Pow_2) + blockDim.x*2);
        return grad + face_grad;
        break;
        }
    default:
        return grad;
        break;
    }
}

/**
 * @brief Compute the total energy of the bond, flatness, bending and dihedral terms from all nodes in the isomer.
 * @param c The constants for the threadIdx^th node.
 * @return Total energy.
*/
INLINE real_t energy(coord3d* X) const {
    BLOCK_SYNC
    real_t arc_energy = (real_t)0.0;

    //(71 + 124) * 3 * N  = 585*N FLOPs
    for (uint8_t j = 0; j < 3; j++ ){
        ArcData arc = ArcData(j, X, node_graph);
        arc_energy += arc.energy(constants);
    }
    switch (T)
    {
    case FLATNESS_ENABLED: {
        FaceData face(X, node_graph);
        return reduction(sdata, arc_energy + face.flatness_energy(constants));
        break;
        }
    case FLAT_BOND: {
        FaceData face(X, node_graph);
        return reduction(sdata, arc_energy + face.flatness_energy(constants));
        break;
        }
    default:
        return reduction(sdata, arc_energy);
        break;
    }
}

INLINE real_t gradnorm(coord3d* X, coord3d& d)const {
    return reduction(sdata, dot(-gradient(X),d));
}

//Bracketing method designed to find upper bound for linesearch method that matches 
//reference python implementation by Buster.
INLINE real_t FindLineSearchBound(coord3d* X, coord3d& r0, coord3d* X1) const{
    real_t bound        = 1e-5;
    bool negative_grad  = true;
    size_t iter         = 0;
    while (negative_grad && iter < 1000)
    {   
        bound *= (real_t)1.5;
        X1[node_id] = X[node_id] + bound * r0;
        real_t gradsum = reduction(sdata, dot(gradient(X1),r0));
        negative_grad = (gradsum < 0);
    }
    return bound;
}


// finds the minimum point of a function f(x) using the bisection method
// X is the current point, X1 is the next point, X2 is the previous point
// r0 is the direction of the line search
// returns the minimum point
INLINE real_t Bisection(coord3d* X, coord3d& r0, coord3d* X1, coord3d* X2){
    real_t dfc = 1; size_t count = 0;
    real_t c; real_t a = 0.0; real_t b = FindLineSearchBound(X,r0,X1);
    coord3d d;
    while (abs(dfc) > 1e-10 && count < 1000){
        count++;
        c =  (a+b)/2;
        X1[node_id] = X[node_id] + c*r0;
        d  =  gradient(X1);
        dfc = reduction(sdata,dot(d,r0)); 

        if (dfc < (real_t)0.0){
            a = c;
        }
        else{
            b = c;
        }
    }
    return c;
}

//Brents Method for line-search using fixed number of iterations.
/**
 * @brief Brent's method for line-search.
 * @param X The coordinates of the nodes.
 * @param r0 The direction of the line-search.
 * @param X1 memory for storing temporary coordinates at a and s.
 * @param X2 memory for storing temporary coordinates at b.
 * @return The step-size s.
*/
INLINE real_t BrentsMethod(coord3d* X, coord3d& r0, coord3d* X1, coord3d* X2)const{
    real_t a,b,s,d;
    a = (real_t)0.0; //b = (real_t)1.0; 

    //To match python reference implementation by Buster.
    b = FindLineSearchBound(X,r0,X1);

    X1[node_id] = X[node_id] + a * r0;
    X2[node_id] = X[node_id] + b * r0;

    real_t f_a = gradnorm(X1,r0);
    real_t f_b = gradnorm(X2,r0);

    if (f_a*f_b > 0)
    {
        return b;
    }
    if (abs(f_a) < abs(f_b)){swap_reals(a,b); swap_reals(f_a,f_b);}

    real_t c = a; real_t f_c = f_a;
    bool flag = true;

    for (uint8_t i = 0; i < 30; i++)
    {   
        // Inverse quadratic interpolation
        if ( (f_a != f_c) && (f_b != f_c) )
        {
            s = a*f_a*f_c / ((f_a - f_b)*(f_a - f_c)) + b * f_a * f_c / ((f_b-f_a)*(f_b-f_c)) + c*f_a*f_b/((f_c-f_a)*(f_c-f_b));
        }else // Secant Method
        {
            s = b - f_b*(b-a)/(f_b-f_a);
        }
        
        bool condition_1 = !(s > (((real_t)3.0*a + b)/(real_t)4.0) && s < b);
        bool condition_2 = flag && (abs(s-b) >= abs(b-c)/(real_t)2.0);
        bool condition_3 = !flag && (abs(s-b) >= abs(c-d)/(real_t)2.0);
        bool condition_4 = flag && (abs(b-c) < (real_t)5e-8);
        bool condition_5 = !flag && (abs(c-d) < (real_t)5e-8);

        if (condition_1 || condition_2 || condition_3 || condition_4 || condition_5)
        {
            s = (a+b) / (real_t)2.0; // Bisection Method
            flag = true;
        }else
        {
            flag = false;
        }
        X1[node_id] = X[node_id] + s * r0;
        real_t f_s = gradnorm(X1,r0);
        d = c;
        c = b; f_c = f_b;
        if (f_a*f_s < 0)
        {
            b = s; f_b = f_s;
        }else
        {
            a = s; f_a = f_s;
        }
        if (abs(f_a) < abs(f_b))
        {
            swap_reals(a,b); swap_reals(f_a,f_b);
        }
    }
    return b;
}

//Golden Section Search, using fixed iterations.
/**
 * @brief Golden Section Search for line-search.
 * @param X The coordinates of the nodes.
 * @param r0 The direction of the line-search.
 * @param X1 memory for storing temporary coordinates at x1.
 * @param X2 memory for storing temporary coordinates at x2.
 * @return The step-size alpha
*/
INLINE real_t GSS(coord3d* X, const coord3d& r0, coord3d* X1, coord3d* X2) const{
    constexpr real_t tau = (real_t)0.6180339887;
    //Line search x - values;
    real_t a = 0.0; real_t b = (real_t)1.0;
    
    real_t x1,  x2;
    x1 = (a + (1 - tau) * (b - a));
    x2 = (a + tau * (b - a));
    //Actual coordinates resulting from each traversal 
    X1[node_id] = X[node_id] + x1 * r0;
    X2[node_id] = X[node_id] + x2 * r0;

    real_t f1 = energy(X1);
    real_t f2 = energy(X2);

    for (uint8_t i = 0; i < 20; i++){
        if (f1 > f2){
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + tau * (b - a);
            X2[node_id] = X[node_id] + x2 * r0;
            f2 = energy(X2);
        }else
        {
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = a + ((real_t)1.0 - tau) * (b - a);
            X1[node_id] = X[node_id] + x1 * r0;
            f1 = energy(X1);
        }
    }
    if (f1 > energy(X)) {return (real_t)0.0;}
    //Line search coefficient
    real_t alpha = (a+b)/(real_t)2.0;
    return alpha;
}

/**
 * @brief Conjugate Gradient Method for energy minimization.
 * @param X The coordinates of the nodes.
 * @param X1 memory for storing temporary coordinates.
 * @param X2 memory for storing temporary coordinates.
 * @param MaxIter The maximum number of iterations.
*/
INLINE  void CG(coord3d* X, coord3d* X1, coord3d* X2, const size_t MaxIter){
    real_t alpha, beta, g0_norm2, s_norm;
    coord3d g0,g1,s;
    g0 = gradient(X);
    s = -g0;
    //Normalize To match reference python implementation by Buster.
    #if USE_MAX_NORM==1
        s_norm = reduction_max(sdata, max(max(s.x,s.y),s.z));
    #else
        s_norm = sqrt(reduction(sdata, dot(s,s)));
    #endif
    s /= s_norm;

    for (size_t i = 0; i < MaxIter; i++){
        alpha = LINESEARCH_METHOD(X,s,X1,X2);
        if (alpha > (real_t)0.0){X1[node_id] = X[node_id] + alpha * s;}
        g1 = gradient(X1);
        //Polak Ribiere method
        g0_norm2 = reduction(sdata, dot(g0, g0));
        beta = d_max(reduction(sdata, dot(g1, (g1 - g0))) / g0_norm2,(real_t)0.0);

        if (alpha > (real_t)0.0){X[node_id] = X1[node_id];}else{ g1 = g0; beta = (real_t) 0.0;}
        s = -g1 + beta*s;
        g0 = g1;
        //Normalize Search Direction using MaxNorm or 2Norm
        #if USE_MAX_NORM==1
            s_norm = reduction_max(sdata, max(max(s.x,s.y),s.z));
        #else
            s_norm = sqrt(reduction(sdata, dot(s,s)));
        #endif
        s /= s_norm;
    }
}
};

/**
 * @brief Checks if the isomer_idx^th isomer has converged, isomer_idx is different for each block.
 * @param B The isomer batch.
 * @param isomer_idx The index of the isomer to check.
 * @param max_iterations The maximum number of iterations, if the isomer has not converged after this many iterations, it is marked as failed.
 * @tparam T The forcefield type.
 * @return void
*/
template <ForcefieldType T>
__device__ void check_batch(IsomerBatch &B, const size_t isomer_idx, const size_t max_iterations){
    DEVICE_TYPEDEFS
    extern __shared__ real_t smem[];
    clear_cache(smem,Block_Size_Pow_2);

    if (isomer_idx < B.isomer_capacity){ //Avoid illegal memory access
    if (B.statuses[isomer_idx] == IsomerStatus::NOT_CONVERGED){
    size_t offset = isomer_idx * blockDim.x;
    Constants constants     = Constants(B, isomer_idx);
    NodeNeighbours node_graph    = NodeNeighbours(B, isomer_idx, smem);
    ForceField FF           = ForceField<T>(node_graph, constants, smem);
    coord3d* X              = reinterpret_cast<coord3d*>(B.X+offset*3);

    coord3d rel_bond_err, rel_angle_err, rel_dihedral_err;
    BLOCK_SYNC
    for (uint8_t j = 0; j < 3; j++){
        auto arc            = ForceField<T>::ArcData(j, X, node_graph);
        #if USE_CONSTANT_INDICES
        d_set(rel_bond_err,      j, abs(abs(arc.bond()       - constants.r0(j))        /constants.r0(j)));
        d_set(rel_angle_err,     j, abs(abs(arc.angle()      - constants.angle0(j))    /constants.angle0(j)));
        d_set(rel_dihedral_err,  j, abs(abs(arc.dihedral()   - constants.inner_dih0(j))/constants.inner_dih0(j)));
        #else
        d_set(rel_bond_err,      j, abs(abs(arc.bond()       - d_get(constants.r0,j))        /d_get(constants.r0,j)));
        d_set(rel_angle_err,     j, abs(abs(arc.angle()      - d_get(constants.angle0,j))    /d_get(constants.angle0,j)));
        d_set(rel_dihedral_err,  j, abs(abs(arc.dihedral()   - d_get(constants.inner_dih0,j))/d_get(constants.inner_dih0,j)));
        #endif
    }

    real_t bond_max         = reduction_max(smem, max(rel_bond_err));
    real_t angle_max        = reduction_max(smem, max(rel_angle_err));
    real_t dihedral_max     = reduction_max(smem, max(rel_dihedral_err));
    real_t bond_rms         = sqrt(reduction(smem,dot(rel_bond_err,rel_bond_err))/blockDim.x);
    real_t angle_rms        = sqrt(reduction(smem,dot(rel_angle_err,rel_angle_err))/blockDim.x);
    real_t dihedral_rms     = sqrt(reduction(smem,dot(rel_dihedral_err,rel_dihedral_err))/blockDim.x);
    real_t bond_mean        = reduction(smem,sum(rel_bond_err))/blockDim.x;
    real_t angle_mean       = reduction(smem,sum(rel_angle_err))/blockDim.x;
    real_t dihedral_mean    = reduction(smem,sum(rel_dihedral_err))/blockDim.x;
    real_t grad_norm        = sqrt(reduction(smem,dot(FF.gradient(X), FF.gradient(X))))/blockDim.x;
    real_t grad_rms         = sqrt(reduction(smem,dot(FF.gradient(X), FF.gradient(X)))/blockDim.x);
    real_t grad_max         = reduction_max(smem, sqrt(dot(FF.gradient(X), FF.gradient(X)))     );
    real_t energy           = FF.energy(X); 
    
    bool converged = ((grad_norm < 1e-2)  && (bond_max < 0.1) && !isnan(grad_norm)) ;
    //if(threadIdx.x + isomer_idx == 0){printf("%d", (int)num_converged); printf("/ %d Fullerenes Converged in Batch \n", (int)gridDim.x);}
    if(threadIdx.x == 0 && B.statuses[isomer_idx] != IsomerStatus::EMPTY){
        if (converged)
        {
            B.statuses[isomer_idx] = IsomerStatus::CONVERGED;
        } else if (B.iterations[isomer_idx] >= max_iterations || isnan(grad_norm)  ) {
            B.statuses[isomer_idx] = IsomerStatus::FAILED;
        }
       
    }}}
}

/**
 * @brief Forcefield Optimizes a batch of isomers.
 * @param B IsomerBatch
 * @param iterations Number of iterations to run
 * @param max_iterations Maximum number of iterations, to compare against, isomers are marked as failed if this is exceeded.
 * @return void
*/
template <ForcefieldType T>
__global__ void optimize_(IsomerBatch B, const size_t iterations, const size_t max_iterations){
    DEVICE_TYPEDEFS
    extern __shared__ real_t smem[];
    clear_cache(smem,Block_Size_Pow_2);
    auto limit = ((B.isomer_capacity + gridDim.x - 1) / gridDim.x ) * gridDim.x;  //Fast ceiling integer division.
    for (int isomer_idx = blockIdx.x; isomer_idx < limit; isomer_idx += gridDim.x){
    BLOCK_SYNC
    if (isomer_idx < B.isomer_capacity){ //Avoid illegal memory access
    if (B.statuses[isomer_idx] == IsomerStatus::NOT_CONVERGED)
    {
        real_t* base_pointer        = smem + Block_Size_Pow_2;
        size_t offset               = isomer_idx * blockDim.x;
        size_t node_id              = threadIdx.x;
        size_t N                    = blockDim.x;


        //Pre-compute force constants and store in registers.
        Constants constants = Constants(B, isomer_idx);
        NodeNeighbours nodeG     = NodeNeighbours(B, isomer_idx, smem);

        //Set VRAM pointer to start of each fullerene, as opposed to at the start of the isomerbatch.
        coord3d* X = reinterpret_cast<coord3d*>(B.X+3*offset);

        //Assign a section of L1 cache to each set of cartesian coordinates X, X1 and X2.
        coord3d* sX =reinterpret_cast<coord3d*>(base_pointer);
        coord3d* X1 =reinterpret_cast<coord3d*>(base_pointer+3*N);
        coord3d* X2 =reinterpret_cast<coord3d*>(base_pointer+6*N);  

                                    
        sX[node_id] = X[node_id];   //Copy cartesian coordinates from DRAM to L1 Cache.
        X           = sX;       //Switch coordinate pointer from DRAM to L1 Cache.



        //Create forcefield struct and use optimization algorithm to optimize the fullerene 
        ForceField FF = ForceField<T>(nodeG, constants, smem);
        FF.CG(X,X1,X2,iterations);
        BLOCK_SYNC
        
        //Copy data back from L1 cache to DRAM 
        reinterpret_cast<coord3d*>(B.X)[offset + threadIdx.x]= X[threadIdx.x];

        if (threadIdx.x == 0) {B.iterations[isomer_idx] += iterations;}
    }}
    //Check the convergence of isomers and assign status accordingly.
    BLOCK_SYNC
    check_batch<T>(B, isomer_idx, max_iterations);
    }
}


#if USE_CONSTANT_INDICES
    #define GET_STAT(fun_1, fun_2, param_fun, equillibrium_param, err_fun) \
            template <ForcefieldType T>\
            __global__ void fun_1(const IsomerBatch B, CuArray<device_real_t> bond_rms){\
                DEVICE_TYPEDEFS\
                extern __shared__ real_t smem[];\
                clear_cache(smem,Block_Size_Pow_2);\
                for (int isomer_idx = blockIdx.x; isomer_idx < B.isomer_capacity; isomer_idx+= gridDim.x){\
                if(B.statuses[isomer_idx] != IsomerStatus::EMPTY){\
                    coord3d rel_err;\
                    size_t offset = isomer_idx * blockDim.x;  \
                    Constants constants     = Constants(B, isomer_idx);\
                    NodeNeighbours node_graph    = NodeNeighbours(B, isomer_idx, smem);\
                    ForceField FF           = ForceField<T>(node_graph, constants, smem);\
                    coord3d* X              = reinterpret_cast<coord3d*>(B.X+offset*3);\
                    for (uint8_t j = 0; j < 3; j++){\
                        auto arc            = ForceField<T>::ArcData(j, X, node_graph);\
                        d_set(rel_err,      j, abs(abs(param_fun       - equillibrium_param(j))        /equillibrium_param(j)));\
                    }\
                    bond_rms.data[isomer_idx]         = err_fun;\
                }}\
            }\
            template <ForcefieldType T>\
            cudaError_t fun_2(const IsomerBatch& B, CuArray<device_real_t>& bond_rms){\
                cudaDeviceSynchronize();\
                cudaSetDevice(B.get_device_id());\
                size_t smem = sizeof(device_coord3d)* (3*B.n_atoms + 4) + sizeof(device_real_t)*Block_Size_Pow_2;\
                static LaunchDims dims((void*)fun_1<T>, B.n_atoms, smem, B.isomer_capacity);\
                dims.update_dims((void*)fun_1<T>, B.n_atoms, smem, B.isomer_capacity);\
                void* kargs[]{(void*)&B, (void*)&bond_rms};\
                auto error = safeCudaKernelCall((void*)fun_1<T>, dims.get_grid(), dims.get_block(), kargs, smem);\
                cudaDeviceSynchronize();\
                return error;\
            }
#else
    #define GET_STAT(fun_1, fun_2, param_fun, equillibrium_param, err_fun) \
        template <ForcefieldType T>\
        __global__ void fun_1(const IsomerBatch B, CuArray<device_real_t> bond_rms){\
            DEVICE_TYPEDEFS\
            extern __shared__ real_t smem[];\
            clear_cache(smem,Block_Size_Pow_2);\
            for (int isomer_idx = blockIdx.x; isomer_idx < B.isomer_capacity; isomer_idx+= gridDim.x){\
            if(B.statuses[isomer_idx] != IsomerStatus::EMPTY){\
                coord3d rel_err;\
                size_t offset = isomer_idx * blockDim.x;  \
                Constants constants     = Constants(B, isomer_idx);\
                NodeNeighbours node_graph    = NodeNeighbours(B, isomer_idx, smem);\
                ForceField FF           = ForceField<T>(node_graph, constants, smem);\
                coord3d* X              = reinterpret_cast<coord3d*>(B.X+offset*3);\
                for (uint8_t j = 0; j < 3; j++){\
                    auto arc            = ForceField<T>::ArcData(j, X, node_graph);\
                    d_set(rel_err,      j, abs(abs(param_fun       - d_get(equillibrium_param,j))        /d_get(equillibrium_param,j)));\
                }\
                bond_rms.data[isomer_idx]         = err_fun;\
            }}\
        }\
        template <ForcefieldType T>\
        cudaError_t fun_2(const IsomerBatch& B, CuArray<device_real_t>& bond_rms){\
            cudaDeviceSynchronize();\
            cudaSetDevice(B.get_device_id());\
            size_t smem = sizeof(device_coord3d)* (3*B.n_atoms + 4) + sizeof(device_real_t)*Block_Size_Pow_2;\
            static LaunchDims dims((void*)fun_1<T>, B.n_atoms, smem, B.isomer_capacity);\
            dims.update_dims((void*)fun_1<T>, B.n_atoms, smem, B.isomer_capacity);\
            void* kargs[]{(void*)&B, (void*)&bond_rms};\
            auto error = safeCudaKernelCall((void*)fun_1<T>, dims.get_grid(), dims.get_block(), kargs, smem);\
            cudaDeviceSynchronize();\
            return error;\
        }
#endif



#define GET_MEAN(fun_1, fun_2, param_fun, err_fun) \
    template <ForcefieldType T>\
    __global__ void fun_1(const IsomerBatch B, CuArray<device_real_t> bond_rms){\
        DEVICE_TYPEDEFS\
        extern __shared__ real_t smem[];\
        clear_cache(smem,Block_Size_Pow_2);\
        for (int isomer_idx = blockIdx.x; isomer_idx < B.isomer_capacity; isomer_idx+= gridDim.x){\
        if(B.statuses[isomer_idx] != IsomerStatus::EMPTY){\
            coord3d rel_err;\
            size_t offset = isomer_idx * blockDim.x;  \
            Constants constants     = Constants(B, isomer_idx);\
            NodeNeighbours node_graph    = NodeNeighbours(B, isomer_idx, smem);\
            ForceField FF           = ForceField<T>(node_graph, constants, smem);\
            coord3d* X              = reinterpret_cast<coord3d*>(B.X+offset*3);\
            for (uint8_t j = 0; j < 3; j++){\
                auto arc            = ForceField<T>::ArcData(j, X, node_graph);\
                d_set(rel_err,      j, abs(param_fun));\
            }\
            bond_rms.data[isomer_idx]         = err_fun;\
        }}\
    }\
    template <ForcefieldType T>\
    cudaError_t fun_2(const IsomerBatch& B, CuArray<device_real_t>& bond_rms){\
        cudaDeviceSynchronize();\
        cudaSetDevice(B.get_device_id());\
        size_t smem = sizeof(device_coord3d)* (3*B.n_atoms + 4) + sizeof(device_real_t)*Block_Size_Pow_2;\
        static LaunchDims dims((void*)fun_1<T>, B.n_atoms, smem, B.isomer_capacity);\
        dims.update_dims((void*)fun_1<T>, B.n_atoms, smem, B.isomer_capacity);\
        void* kargs[]{(void*)&B, (void*)&bond_rms};\
        auto error = safeCudaKernelCall((void*)fun_1<T>, dims.get_grid(), dims.get_block(), kargs, smem);\
        cudaDeviceSynchronize();\
        return error;\
    }

#define GET_RRMSE(fun_1, fun_2, param_fun, err_fun) \
    template <ForcefieldType T>\
    __global__ void fun_1(const IsomerBatch B, CuArray<device_real_t> bond_rms){\
        DEVICE_TYPEDEFS\
        extern __shared__ real_t smem[];\
        clear_cache(smem,Block_Size_Pow_2);\
        for (int isomer_idx = blockIdx.x; isomer_idx < B.isomer_capacity; isomer_idx+= gridDim.x){\
        if(B.statuses[isomer_idx] != IsomerStatus::EMPTY){\
            coord3d top;\
            coord3d bot;\
            size_t offset = isomer_idx * blockDim.x;  \
            Constants constants     = Constants(B, isomer_idx);\
            NodeNeighbours node_graph    = NodeNeighbours(B, isomer_idx, smem);\
            ForceField FF           = ForceField<T>(node_graph, constants, smem);\
            coord3d* X              = reinterpret_cast<coord3d*>(B.X+offset*3);\
            for (uint8_t j = 0; j < 3; j++){\
                auto arc            = ForceField<T>::ArcData(j, X, node_graph);\
                d_set(top,      j, param_fun -  err_fun);\
                d_set(bot,      j, err_fun);\
            }\
            bond_rms.data[isomer_idx]         = sqrt((reduction(smem, dot(top,top))/reduction(smem,dot(bot,bot)))/(blockDim.x*3));\
        }}\
    }\
    template <ForcefieldType T>\
    cudaError_t fun_2(const IsomerBatch& B, CuArray<device_real_t>& bond_rms){\
        cudaDeviceSynchronize();\
        cudaSetDevice(B.get_device_id());\
        size_t smem = sizeof(device_coord3d)* (3*B.n_atoms + 4) + sizeof(device_real_t)*Block_Size_Pow_2;\
        static LaunchDims dims((void*)fun_1<T>, B.n_atoms, smem, B.isomer_capacity);\
        dims.update_dims((void*)fun_1<T>, B.n_atoms, smem, B.isomer_capacity);\
        void* kargs[]{(void*)&B, (void*)&bond_rms};\
        auto error = safeCudaKernelCall((void*)fun_1<T>, dims.get_grid(), dims.get_block(), kargs, smem);\
        cudaDeviceSynchronize();\
        return error;\
    }

#define GET_RMSE(fun_1, fun_2, param_fun, err_fun) \
    template <ForcefieldType T>\
    __global__ void fun_1(const IsomerBatch B, CuArray<device_real_t> bond_rms){\
        DEVICE_TYPEDEFS\
        extern __shared__ real_t smem[];\
        clear_cache(smem,Block_Size_Pow_2);\
        for (int isomer_idx = blockIdx.x; isomer_idx < B.isomer_capacity; isomer_idx+= gridDim.x){\
        if(B.statuses[isomer_idx] != IsomerStatus::EMPTY){\
            coord3d top;\
            size_t offset = isomer_idx * blockDim.x;  \
            Constants constants     = Constants(B, isomer_idx);\
            NodeNeighbours node_graph    = NodeNeighbours(B, isomer_idx, smem);\
            ForceField FF           = ForceField<T>(node_graph, constants, smem);\
            coord3d* X              = reinterpret_cast<coord3d*>(B.X+offset*3);\
            for (uint8_t j = 0; j < 3; j++){\
                auto arc            = ForceField<T>::ArcData(j, X, node_graph);\
                d_set(top,      j, param_fun - err_fun);\
            }\
            bond_rms.data[isomer_idx]         = sqrt(reduction(smem, dot(top,top))/(blockDim.x*3));\
        }}\
    }\
    template <ForcefieldType T>\
    cudaError_t fun_2(const IsomerBatch& B, CuArray<device_real_t>& bond_rms){\
        cudaDeviceSynchronize();\
        cudaSetDevice(B.get_device_id());\
        size_t smem = sizeof(device_coord3d)* (3*B.n_atoms + 4) + sizeof(device_real_t)*Block_Size_Pow_2;\
        static LaunchDims dims((void*)fun_1<T>, B.n_atoms, smem, B.isomer_capacity);\
        dims.update_dims((void*)fun_1<T>, B.n_atoms, smem, B.isomer_capacity);\
        void* kargs[]{(void*)&B, (void*)&bond_rms};\
        auto error = safeCudaKernelCall((void*)fun_1<T>, dims.get_grid(), dims.get_block(), kargs, smem);\
        cudaDeviceSynchronize();\
        return error;\
    }

#define GET_INTERNAL(fun_1, fun_2, param_fun) \
    __global__ void fun_1(const IsomerBatch B, CuArray<device_real_t> bond_rms){\
        DEVICE_TYPEDEFS\
        extern __shared__ real_t smem[];\
        clear_cache(smem,Block_Size_Pow_2);\
        for (int isomer_idx = blockIdx.x; isomer_idx < B.isomer_capacity; isomer_idx+= gridDim.x){\
        if(B.statuses[isomer_idx] != IsomerStatus::EMPTY){\
            coord3d rel_err;\
            size_t offset = isomer_idx * blockDim.x;  \
            Constants constants     = Constants(B, isomer_idx);\
            NodeNeighbours node_graph    = NodeNeighbours(B, isomer_idx, smem);\
            ForceField FF           = ForceField<FORCEFIELD_VERSION>(node_graph, constants, smem);\
            coord3d* X              = reinterpret_cast<coord3d*>(B.X+offset*3);\
            for (uint8_t j = 0; j < 3; j++){\
                auto arc            = ForceField<FORCEFIELD_VERSION>::ArcData(j, X, node_graph);\
                d_set(rel_err,      j, param_fun);\
            }\
            reinterpret_cast<coord3d*>(bond_rms.data)[isomer_idx*B.n_atoms + threadIdx.x]         =  {rel_err.x, rel_err.y, rel_err.z};\
        }}\
    }\
    cudaError_t fun_2(const IsomerBatch& B, CuArray<device_real_t>& bond_rms){\
        cudaDeviceSynchronize();\
        cudaSetDevice(B.get_device_id());\
        size_t smem = sizeof(device_coord3d)* (3*B.n_atoms + 4) + sizeof(device_real_t)*Block_Size_Pow_2;\
        static LaunchDims dims((void*)fun_1, B.n_atoms, smem, B.isomer_capacity);\
        dims.update_dims((void*)fun_1, B.n_atoms, smem, B.isomer_capacity);\
        void* kargs[]{(void*)&B, (void*)&bond_rms};\
        auto error = safeCudaKernelCall((void*)fun_1, dims.get_grid(), dims.get_block(), kargs, smem);\
        cudaDeviceSynchronize();\
        return error;\
    }

GET_INTERNAL(get_bonds_,get_bonds, arc.bond())
GET_INTERNAL(get_angles_,get_angles, arc.angle())
GET_INTERNAL(get_dihedrals_,get_dihedrals, arc.dihedral())

GET_RRMSE(get_bond_rrmse_,get_bond_rrmse, arc.bond(), d_get(constants.r0,j))
GET_RRMSE(get_angle_rrmse_,get_angle_rrmse, acos(arc.angle()), acos(d_get(constants.angle0,j)))
GET_RRMSE(get_dihedral_rrmse_,get_dihedral_rrmse, acos(arc.dihedral()), acos(d_get(constants.inner_dih0,j)))

GET_RMSE(get_bond_rmse_,get_bond_rmse, arc.bond(), d_get(constants.r0,j))
GET_RMSE(get_angle_rmse_,get_angle_rmse, acos(arc.angle()), acos(d_get(constants.angle0,j)))
GET_RMSE(get_dihedral_rmse_,get_dihedral_rmse, acos(arc.dihedral()), acos(d_get(constants.inner_dih0,j)))

GET_STAT(get_bond_max_,get_bond_max, arc.bond(), constants.r0, reduction_max(smem, max(rel_err)))
GET_STAT(get_angle_rms_,get_angle_rms, arc.angle(), constants.angle0, sqrt(reduction(smem,dot(rel_err,rel_err))/blockDim.x);)
GET_STAT(get_angle_max_,get_angle_max, arc.angle(), constants.angle0, reduction_max(smem, max(rel_err)))
GET_STAT(get_dihedral_max_,get_dihedral_max, arc.dihedral(), constants.inner_dih0, reduction_max(smem, max(rel_err)))
GET_STAT(get_energies_,get_energies, arc.dihedral(), constants.inner_dih0, FF.energy(X))
GET_STAT(get_gradient_norm_,get_gradient_norm, arc.dihedral(), constants.inner_dih0, sqrt(reduction(smem,dot(FF.gradient(X), FF.gradient(X))))/blockDim.x)
GET_STAT(get_gradient_rms_,get_gradient_rms, arc.dihedral(), constants.inner_dih0, sqrt(reduction(smem,dot(FF.gradient(X), FF.gradient(X)))/blockDim.x))
GET_STAT(get_gradient_max_,get_gradient_max, arc.dihedral(), constants.inner_dih0, reduction_max(smem, sqrt(dot(FF.gradient(X), FF.gradient(X)))     ))

GET_MEAN(get_bond_mae_, get_bond_mae, abs(arc.bond() - d_get(constants.r0,j)), reduction(smem,sum(rel_err)/3.0f)/blockDim.x)
GET_MEAN(get_angle_mae_, get_angle_mae, abs(acos(arc.angle()) - acos(d_get(constants.angle0,j))), reduction(smem,sum(rel_err)/3.0f)/blockDim.x)
GET_MEAN(get_dihedral_mae_, get_dihedral_mae, abs(acos(arc.dihedral()) - acos(d_get(constants.inner_dih0,j))), reduction(smem,sum(rel_err)/3.0f)/blockDim.x)

GET_MEAN(get_energy_, get_energy, arc.bond(), FF.energy(X)/blockDim.x)
GET_MEAN(get_bond_mean_,get_bond_mean, arc.bond(), reduction(smem,sum(rel_err)/3.0f)/blockDim.x)
GET_MEAN(get_angle_mean_,get_angle_mean, arc.angle(), reduction(smem,sum(rel_err)/3.0f)/blockDim.x)
GET_MEAN(get_dihedral_mean_,get_dihedral_mean, arc.dihedral(), reduction(smem,sum(rel_err)/3.0f)/blockDim.x)
GET_MEAN(get_gradient_mean_,get_gradient_mean, d_get(FF.gradient(X),j), reduction(smem,sum(rel_err)/3.0f)/blockDim.x)
GET_MEAN(get_flat_mean_,get_flat_mean, d_get(FF.gradient(X),j), reduction(smem,ForceField<FORCEFIELD_VERSION>::FaceData(X, node_graph).flatness())/(blockDim.x/2 + 2 ))
GET_MEAN(get_flat_rmse_,get_flat_rmse, d_get(FF.gradient(X),j), sqrt(reduction(smem,ForceField<FORCEFIELD_VERSION>::FaceData(X, node_graph).flatness() * ForceField<FORCEFIELD_VERSION>::FaceData(X, node_graph).flatness())/(blockDim.x/2 + 2 )) )
GET_MEAN(get_flat_max_,get_flat_max, d_get(FF.gradient(X),j), reduction_max(smem,ForceField<FORCEFIELD_VERSION>::FaceData(X, node_graph).flatness())   )

int optimal_batch_size(const int N, const int device_id) {
    cudaSetDevice(device_id);
    static size_t smem = sizeof(device_coord3d)*3*N + sizeof(device_real_t)*Block_Size_Pow_2;
    static LaunchDims dims((void*)optimize_<FORCEFIELD_VERSION>, N, smem);
    dims.update_dims((void*)optimize_<FORCEFIELD_VERSION>, N, smem);
    return dims.get_grid().x;
}

float kernel_time = 0.0;
std::chrono::microseconds time_spent(){
    return std::chrono::microseconds((int) (kernel_time*1000.f));
}


template <ForcefieldType T>
cudaError_t optimize(IsomerBatch& B, const size_t iterations, const size_t max_iterations, const LaunchCtx& ctx, const LaunchPolicy policy){
    cudaSetDevice(B.get_device_id());
    static std::vector<bool> first_call(16, true);
    static cudaEvent_t start[16], stop[16];
    float single_kernel_time = 0.0;
    auto dev = B.get_device_id();
    if(first_call[dev]) {cudaEventCreate(&start[dev]); cudaEventCreate(&stop[dev]);}
        
    //If launch ploicy is synchronous then wait.
    if(policy == LaunchPolicy::SYNC) {ctx.wait();}
    else if(policy == LaunchPolicy::ASYNC && !first_call[dev]){
        //Records time from previous kernel call
        cudaEventElapsedTime(&single_kernel_time, start[dev], stop[dev]);
        kernel_time += single_kernel_time;
    }

    size_t smem = sizeof(device_coord3d)* (3*B.n_atoms + 4) + sizeof(device_real_t)*Block_Size_Pow_2;
    static LaunchDims dims((void*)optimize_<T>, B.n_atoms, smem, B.isomer_capacity);
    dims.update_dims((void*)optimize_<T>, B.n_atoms, smem, B.isomer_capacity);
    void* kargs[]{(void*)&B, (void*)&iterations, (void*)&max_iterations};

    cudaEventRecord(start[dev], ctx.stream);
    auto error = safeCudaKernelCall((void*)optimize_<T>, dims.get_grid(), dims.get_block(), kargs, smem, ctx.stream);
    cudaEventRecord(stop[dev], ctx.stream);
    
    if(policy == LaunchPolicy::SYNC) {
        ctx.wait();
        cudaEventElapsedTime(&single_kernel_time, start[dev], stop[dev]);
        kernel_time += single_kernel_time;
    }
    printLastCudaError("Forcefield: ");
    first_call[dev] = false;
    return error;
}

int declare_generics(){
    IsomerBatch B(20,1,DEVICE_BUFFER);
    CuArray<device_real_t> arr(1);

    optimize<BUSTER>(B,100,100);
    get_angle_max<BUSTER>(B,arr);
    get_bond_max<BUSTER>(B,arr);
    get_dihedral_max<BUSTER>(B,arr);
    get_angle_mae<BUSTER>(B,arr);
    get_bond_mae<BUSTER>(B,arr);
    get_dihedral_mae<BUSTER>(B,arr);
    get_angle_rrmse<BUSTER>(B,arr);
    get_bond_rrmse<BUSTER>(B,arr);
    get_dihedral_rrmse<BUSTER>(B,arr);
    get_angle_rmse<BUSTER>(B,arr);
    get_bond_rmse<BUSTER>(B,arr);
    get_dihedral_rmse<BUSTER>(B,arr);
    get_angle_mean<BUSTER>(B,arr);
    get_bond_mean<BUSTER>(B,arr);
    get_flat_mean<BUSTER>(B,arr);
    get_flat_max<BUSTER>(B,arr);
    get_flat_rmse<BUSTER>(B,arr);
    get_dihedral_mean<BUSTER>(B,arr);
    get_gradient_max<BUSTER>(B,arr);
    get_gradient_rms<BUSTER>(B,arr);
    get_gradient_mean<BUSTER>(B,arr);
    get_gradient_norm<BUSTER>(B,arr);
    get_energies<BUSTER>(B,arr);

    optimize<FLATNESS_ENABLED>(B,100,100);
    get_angle_max<FLATNESS_ENABLED>(B,arr);
    get_bond_max<FLATNESS_ENABLED>(B,arr);
    get_dihedral_max<FLATNESS_ENABLED>(B,arr);
    get_angle_mae<FLATNESS_ENABLED>(B,arr);
    get_bond_mae<FLATNESS_ENABLED>(B,arr);
    get_dihedral_mae<FLATNESS_ENABLED>(B,arr);
    get_angle_rrmse<FLATNESS_ENABLED>(B,arr);
    get_bond_rrmse<FLATNESS_ENABLED>(B,arr);
    get_dihedral_rrmse<FLATNESS_ENABLED>(B,arr);
    get_angle_rmse<FLATNESS_ENABLED>(B,arr);
    get_bond_rmse<FLATNESS_ENABLED>(B,arr);
    get_dihedral_rmse<FLATNESS_ENABLED>(B,arr);
    get_angle_mean<FLATNESS_ENABLED>(B,arr);
    get_bond_mean<FLATNESS_ENABLED>(B,arr);
    get_dihedral_mean<FLATNESS_ENABLED>(B,arr);
    get_flat_mean<FLATNESS_ENABLED>(B,arr);
    get_flat_max<FLATNESS_ENABLED>(B,arr);
    get_flat_rmse<FLATNESS_ENABLED>(B,arr);
    get_gradient_max<FLATNESS_ENABLED>(B,arr);
    get_gradient_rms<FLATNESS_ENABLED>(B,arr);
    get_gradient_mean<FLATNESS_ENABLED>(B,arr);
    get_gradient_norm<FLATNESS_ENABLED>(B,arr);
    get_energies<FLATNESS_ENABLED>(B,arr);

    optimize<WIRZ>(B,100,100);
    get_angle_max<WIRZ>(B,arr);
    get_bond_max<WIRZ>(B,arr);
    get_dihedral_max<WIRZ>(B,arr);
    get_angle_mae<WIRZ>(B,arr);
    get_bond_mae<WIRZ>(B,arr);
    get_dihedral_mae<WIRZ>(B,arr);
    get_angle_rrmse<WIRZ>(B,arr);
    get_bond_rrmse<WIRZ>(B,arr);
    get_dihedral_rrmse<WIRZ>(B,arr);
    get_angle_rmse<WIRZ>(B,arr);
    get_bond_rmse<WIRZ>(B,arr);
    get_dihedral_rmse<WIRZ>(B,arr);
    get_angle_mean<WIRZ>(B,arr);
    get_bond_mean<WIRZ>(B,arr);
    get_dihedral_mean<WIRZ>(B,arr);
    get_flat_mean<WIRZ>(B,arr);
    get_flat_max<WIRZ>(B,arr);
    get_flat_rmse<WIRZ>(B,arr);
    get_gradient_max<WIRZ>(B,arr);
    get_gradient_rms<WIRZ>(B,arr);
    get_gradient_mean<WIRZ>(B,arr);
    get_gradient_norm<WIRZ>(B,arr);
    get_energies<WIRZ>(B,arr);

    optimize<FLAT_BOND>(B,100,100);
    get_angle_max<FLAT_BOND>(B,arr);
    get_bond_max<FLAT_BOND>(B,arr);
    get_dihedral_max<FLAT_BOND>(B,arr);
    get_angle_mae<FLAT_BOND>(B,arr);
    get_bond_mae<FLAT_BOND>(B,arr);
    get_dihedral_mae<FLAT_BOND>(B,arr);
    get_angle_rrmse<FLAT_BOND>(B,arr);
    get_bond_rrmse<FLAT_BOND>(B,arr);
    get_dihedral_rrmse<FLAT_BOND>(B,arr);
    get_angle_rmse<FLAT_BOND>(B,arr);
    get_bond_rmse<FLAT_BOND>(B,arr);
    get_dihedral_rmse<FLAT_BOND>(B,arr);
    get_angle_mean<FLAT_BOND>(B,arr);
    get_bond_mean<FLAT_BOND>(B,arr);
    get_dihedral_mean<FLAT_BOND>(B,arr);
    get_flat_mean<FLAT_BOND>(B,arr);
    get_flat_max<FLAT_BOND>(B,arr);
    get_flat_rmse<FLAT_BOND>(B,arr);
    get_gradient_max<FLAT_BOND>(B,arr);
    get_gradient_rms<FLAT_BOND>(B,arr);
    get_gradient_mean<FLAT_BOND>(B,arr);
    get_gradient_norm<FLAT_BOND>(B,arr);
    get_energies<FLAT_BOND>(B,arr);

    return 1;
}


}}
