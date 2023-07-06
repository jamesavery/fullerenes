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
namespace isomerspace_hessian{
#include "device_includes.cu"

// This struct was made to reduce signature cluttering of device functions, it is simply a container for default arguments which are shared between functions
template <ForcefieldType T>
struct ForceField{
    DEVICE_TYPEDEFS;
    
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
            auto a = Xf[0][0] * Xf[0][0] + Xf[1][0] * Xf[1][0] + Xf[2][0] * Xf[2][0] + Xf[3][0] * Xf[3][0] + Xf[4][0] * Xf[4][0] + Xf[5][0] * Xf[5][0],
                 b = Xf[0][0] * Xf[0][1] + Xf[1][0] * Xf[1][1] + Xf[2][0] * Xf[2][1] + Xf[3][0] * Xf[3][1] + Xf[4][0] * Xf[4][1] + Xf[5][0] * Xf[5][1],
                 c = Xf[0][0] * Xf[0][2] + Xf[1][0] * Xf[1][2] + Xf[2][0] * Xf[2][2] + Xf[3][0] * Xf[3][2] + Xf[4][0] * Xf[4][2] + Xf[5][0] * Xf[5][2],
                 d = Xf[0][1] * Xf[0][1] + Xf[1][1] * Xf[1][1] + Xf[2][1] * Xf[2][1] + Xf[3][1] * Xf[3][1] + Xf[4][1] * Xf[4][1] + Xf[5][1] * Xf[5][1],
                 e = Xf[0][1] * Xf[0][2] + Xf[1][1] * Xf[1][2] + Xf[2][1] * Xf[2][2] + Xf[3][1] * Xf[3][2] + Xf[4][1] * Xf[4][2] + Xf[5][1] * Xf[5][2],
                 f = Xf[0][2] * Xf[0][2] + Xf[1][2] * Xf[1][2] + Xf[2][2] * Xf[2][2] + Xf[3][2] * Xf[3][2] + Xf[4][2] * Xf[4][2] + Xf[5][2] * Xf[5][2];
            //Xf * Xf^T In closed form.
            A = symMat3(a,b,c,d,e,f);

            //A is positive-semi-definite so all eigenvalues are non-negative
            lambdas = A.eigenvalues();
            lambda_f = d_min(d_min(lambdas[0], lambdas[1]), lambdas[2]);  
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
        return c.f_flat * harmonic_energy(flatness(),(real_t)0.);
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
	  auto lam = lambda_f;
	  centroids[threadIdx.x] = centroid;
	  norms[threadIdx.x] = A.eigenvector3x3(lam);
        }
        BLOCK_SYNC

        coord3d grad = {(real_t)0., (real_t)0., (real_t)0.};
        for(unsigned char j = 0; j < 3; j++) grad += dot(Xa - centroids[d_get(face_neighbours,j)], norms[d_get(face_neighbours,j)]) * norms[d_get(face_neighbours,j)];
        return c.f_flat * (real_t)2. * grad;
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
        coord3d ap, am, ab, ac, ad, mp, db, bc;
        coord3d X_a = X[a], X_b = X[d_get(G.cubic_neighbours,j)];
        //Compute the arcs ab, ac, ad, bp, bm, ap, am, mp, bc and cd
        ab = (X_b - X_a);  rabn = bond_length(ab); abh = rabn * ab;
        ac = (X[d_get(G.cubic_neighbours,(j+1)%3)] - X_a); racn = bond_length(ac); ach = racn * ac; rab = non_resciprocal_bond_length(ab);
        ad = (X[d_get(G.cubic_neighbours,(j+2)%3)] - X_a); radn = bond_length(ad); adh = radn * ad;
        db = (X_b - X[d_get(G.cubic_neighbours,(j+2)%3)]); rdbn = bond_length(db); dbh = rdbn * db;
        coord3d bp = (X[d_get(G.next_on_face,j)] - X_b); 
        coord3d bm = (X[d_get(G.prev_on_face,j)] - X_b); 
        
        rbpn = bond_length(bp); bph = bp * rbpn;
        rbmn = bond_length(bm); bmh = bm * rbmn;

        ap = bp + ab; rapn = bond_length(ap); aph = rapn * ap;
        am = bm + ab; ramn = bond_length(am); amh = ramn * am;
        mp = bp - bm; rmpn = bond_length(mp); mph = rmpn * mp;
        bc = ac - ab; rbcn = bond_length(bc); bch = rbcn * bc;
        cdh = unit_vector(ad - ac);

        //Compute inverses of some arcs, these are subject to be omitted if the equations are adapted appropriately with inversion of signs.
        bah = -abh;
        mbh = -bmh;
        pah = -aph;
        pbh = -bph;
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

    INLINE mat3 harmonic_energy_hessian(const real_t p0, const real_t p, const coord3d grad_a, const coord3d grad_b, const mat3& hessp) const{
        return hessp*(p-p0) + tensor_product(grad_a,grad_b);
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
    INLINE real_t angle() const {return dot(abh,ach);}

    INLINE real_t normalized_angle_err() const {return acos((float)(float)dot(abh,ach));}

    //Returns outer angle m, used only diagnostically.
    INLINE real_t outer_angle_m() const {return -dot(abh, bmh);} //Compute outer angle. ab,bm

    //Returns outer angle p, used only diagnostically.
    INLINE real_t outer_angle_p() const{return -dot(abh, bph);} //Compute outer angle. ab,bp

    //Returns the inner dihedral angle for the current arc. Used here only for energy calculation, 
    //otherwise embedded in dihedral computation because the planes and angles that make up the dihedral angle computation are required for derivative computation.
    //50 FLOPs
    /**
     * @brief Compute the dihedral angle between the planes (abc,bcd), (acd,bcd) and (abd,bcd). For j = 0, 1 or 2 respectively. 
     * @return The dihedral angle.
    */
    INLINE real_t dihedral() const 
    { 
        coord3d nabc, nbcd; real_t cos_b, cos_c, r_sin_b, r_sin_c;
        cos_b = dot(bah,bch);  r_sin_b = (real_t)1.0/SQRT((real_t)1.0 - cos_b*cos_b); nabc = cross(bah, bch) * r_sin_b;
        cos_c = dot(-bch,cdh); r_sin_c = (real_t)1.0/SQRT((real_t)1.0 - cos_c*cos_c); nbcd = cross(-bch,cdh) * r_sin_c;
        return dot(nabc, nbcd);
    }
    //Returns the Outer-dihedral-a wrt. current arc, only accessed diagnostically (internal coordinate).
    /**
     * @brief Compute the dihedral angle between the planes $(b-a-b_m, a-b_m-b_p)$, $(c-a-c_m, a-c_m-c_p)$ and $(d-a-d_m, a-d_m-d_p)$. For j = 0, 1 or 2 respectively.
     * @return The dihedral angle.
    */
    INLINE real_t outer_dihedral_a() const
    {
        coord3d nbam_hat, namp_hat; real_t cos_a, cos_m, r_sin_a, r_sin_m;
        cos_a = dot(abh,amh); r_sin_a = (real_t)1.0/SQRT((real_t)1.0 - cos_a*cos_a); nbam_hat = cross(abh,amh) * r_sin_a;
        cos_m = dot(-amh,mph); r_sin_m = (real_t)1.0/SQRT((real_t)1.0 - cos_m*cos_m); namp_hat = cross(-amh,mph) * r_sin_m;
        real_t cos_beta = dot(nbam_hat, namp_hat); //Outer Dihedral angle bam, amp
        return cos_beta;
    }
    //Returns the Outer-dihedral-m wrt. current arc, only accessed diagnostically (internal coordinate).
    /**
     * @brief Compute the dihedral angle between the planes $(b-b_m-b_p, b_m-b_p-a)$, $(c-c_m-c_p, c_m-c_p-a)$ and $(d-d_m-d_p, d_m-d_p-a)$. For j = 0, 1 or 2 respectively.
     * @return The dihedral angle.
    */
    INLINE real_t outer_dihedral_m() const
    {
        coord3d nbmp_hat, nmpa_hat; real_t cos_m, cos_p, r_sin_m, r_sin_p;
        cos_m = dot(mbh,mph);  r_sin_m = (real_t)1.0/SQRT((real_t)1.0 - cos_m*cos_m); nbmp_hat = cross(mbh,mph) * r_sin_m;
        cos_p = dot(-mph,pah); r_sin_p = (real_t)1.0/SQRT((real_t)1.0 - cos_p*cos_p); nmpa_hat = cross(-mph,pah) * r_sin_p;
        //Cosine to the outer dihedral angle constituted by the planes bmp and mpa
        real_t cos_beta = dot(nbmp_hat, nmpa_hat); //Outer dihedral angle bmp,mpa.
        return cos_beta;    
    }
    //Returns the Outer-dihedral-p wrt. current arc, only accessed diagnostically (internal coordinate).
    /**
     * @brief Compute the dihedral angle between the planes $(b-b_p-a, b_p-a-b_m)$, $(c-c_p-a, c_p-a-c_m)$ and $(d-d_p-a, d_p-a-d_m)$. For j = 0, 1 or 2 respectively.
     * @return The dihedral angle.
    */
    INLINE real_t outer_dihedral_p() const
    {
        coord3d nbpa_hat, npam_hat; real_t cos_p, cos_a, r_sin_p, r_sin_a;
        cos_a = dot(aph,amh);  r_sin_a = (real_t)1.0/SQRT((real_t)1.0 - cos_a*cos_a); npam_hat = cross(aph,amh)  * r_sin_a;
        cos_p = dot(pbh,-aph); r_sin_p = (real_t)1.0/SQRT((real_t)1.0 - cos_p*cos_p); nbpa_hat = cross(pbh,-aph) * r_sin_p;
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
        coord3d grad = cos_angle * (abh * rabn + ach * racn) - abh * racn - ach* rabn; //Derivative of inner angle: Eq. 21. 
        return d_get(c.f_inner_angle,j) * harmonic_energy_gradient(d_get(c.angle0,j), cos_angle, grad); //Harmonic Energy Gradient: Eq. 21. multiplied by harmonic term.
    }
    
    INLINE mat3 bond_hessian_a(const Constants& c) const{
        auto grad_a = -abh;
        auto hessp = (identity3() - tensor_product(abh,abh))*rabn;
        return d_get(c.f_bond,j) * harmonic_energy_hessian(d_get(c.r0,j), rab, grad_a, grad_a, hessp);   
    }

    INLINE mat3 bond_hessian_b(const Constants& c) const{
        auto grad_a = -abh;
        auto grad_b = abh;
        auto hessp = (tensor_product(abh,abh) - identity3())*rabn;
        return d_get(c.f_bond,j) * harmonic_energy_hessian(d_get(c.r0,j), rab, grad_a, grad_b, hessp);
    }

    INLINE mat3 inner_angle_hessian_a(const Constants& c) const{
        real_t cos_angle = angle(); //Inner angle of arcs ab,ac.
        coord3d grad_a = (abh * angle() - ach) * rabn + (ach * angle() - abh) * racn; 
        
        //TensorProduct[abh, 1/abn * (abh*cost - ach) + 1/acn * (ach*cost-abh)] + cost/abn * (TensorProduct[abh,abh] - IdentityMatrix[3]) - 1/acn * (TensorProduct[ach,ach] - IdentityMatrix[3])
        auto G = tensor_product(abh, rabn * (abh * angle() - ach) + racn * (ach * angle() - abh))
                + angle()*rabn * (tensor_product(abh, abh) - identity3()) 
                - racn * (tensor_product(ach, ach) - identity3());

        //F := TensorProduct[ach, 1/abn * (abh*cost - ach) + 1/acn * (ach*cost-abh)] + cost/acn * (TensorProduct[ach,ach] - IdentityMatrix[3]) - 1/abn * (TensorProduct[abh,abh] - IdentityMatrix[3])
        auto F = tensor_product(ach, rabn * (abh * angle() - ach) + racn * (ach * angle() - abh)) + angle()*racn * (tensor_product(ach, ach) - identity3()) - rabn * (tensor_product(abh, abh) - identity3());

        auto P1 = tensor_product(abh * angle() - ach, abh * rabn*rabn);
        auto P2 = rabn*G;
        auto P3 = tensor_product(ach * angle() - abh, ach * racn*racn);
        auto P4 = racn*F;
        return d_get(c.f_inner_angle,j) * harmonic_energy_hessian(d_get(c.angle0,j), cos_angle, grad_a, grad_a, P1+P2+P3+P4); //Harmonic Energy Hessian
    }

    INLINE mat3 inner_angle_hessian_b(const Constants& c) const{
        coord3d grad_a = (abh * angle() - ach) * rabn + (ach * angle() - abh) * racn; 
        coord3d grad_b = (ach - abh * angle()) * rabn;
        auto G = tensor_product(abh, rabn * (ach - abh * angle())) + angle()*rabn * (identity3() - tensor_product(abh, abh));
        auto F = tensor_product(ach, rabn * (ach - abh * angle())) - rabn * (identity3() - tensor_product(abh, abh));
        auto P1 = tensor_product(abh * angle() - ach, -abh * rabn*rabn);
        
        auto P2 = rabn*G;
        auto P4 = racn*F;
        return d_get(c.f_inner_angle,j) * harmonic_energy_hessian(d_get(c.angle0,j), angle(), grad_a, grad_b, P1+P2+P4); //Harmonic Energy Hessian
    }

    INLINE mat3 inner_angle_hessian_c(const Constants& c) const{
        auto grad_a = (abh * angle() - ach) * rabn + (ach * angle() - abh) * racn;
        auto grad_c = (abh - angle() * ach) * racn;
        auto G = tensor_product(abh, racn * (abh - ach * angle())) - racn * (identity3() - tensor_product(ach, ach));
        auto F = tensor_product(ach, racn * (abh - ach * angle())) + angle()*racn * (identity3() - tensor_product(ach, ach));
        auto P2 = rabn*G;
        auto P3 = tensor_product(ach * angle() - abh, -ach * racn*racn);
        auto P4 = racn*F;
        return d_get(c.f_inner_angle,j) * harmonic_energy_hessian(d_get(c.angle0,j), angle(), grad_a, grad_c, P2+P3+P4); //Harmonic Energy Hessian
    }

    INLINE mat3 outer_angle_hessian_m_a(const Constants& c) const{
        auto cost = dot(bah, bmh); //Compute outer angle. ab,bm
        auto gradba = rabn * (identity3() - tensor_product(bah, bah));
        auto grad_a = (bmh - bah * cost) * rabn;
        auto P1 = tensor_product(bmh - bah * cost, -bah * rabn*rabn);
        auto P2 = -rabn * (tensor_product(bah, grad_a) + cost * gradba);
        return d_get(c.f_outer_angle_m,j) * harmonic_energy_hessian(d_get(c.outer_angle_m0,j), cost, grad_a, grad_a, P1+P2); //Harmonic Energy Hessian
    }

    INLINE mat3 outer_angle_hessian_m_b(const Constants& c) const{
        auto cost = dot(bah, bmh); //Compute outer angle. ba,bm
        auto gradba = rabn * (tensor_product(bah, bah) - identity3());
        auto gradbm = rbmn * (tensor_product(bmh, bmh) - identity3());
        auto grad_b = rabn * (bah * cost - bmh) + rbmn * (bmh * cost - bah);
        auto grad_a = (bmh - bah * cost) * rabn;
        auto P1 = tensor_product(bmh - bah * cost, bah * rabn*rabn);
        auto P3 = rabn * (gradbm - (tensor_product(bah, grad_b) + cost * gradba));
        return d_get(c.f_outer_angle_m,j) * harmonic_energy_hessian(d_get(c.outer_angle_m0,j), cost, grad_a, grad_b, P1+P3); //Harmonic Energy Hessian
    }

    INLINE mat3 outer_angle_hessian_m_m(const Constants& c) const{
        auto cost = dot(bah, bmh); //Compute outer angle. ba,bm
        auto gradbm = rbmn * (identity3() - tensor_product(bmh, bmh));
        auto grad_a = (bmh - bah * cost) * rabn;
        auto grad_m = rbmn * (bah - bmh * cost);
        auto P1 = rabn * (gradbm - tensor_product(bah, grad_m));
        return d_get(c.f_outer_angle_m,j) * harmonic_energy_hessian(d_get(c.outer_angle_m0,j), cost, grad_a, grad_m, P1); //Harmonic Energy Hessian
    }

    INLINE mat3 outer_angle_hessian_p_a(const Constants& c) const{
        auto cost = dot(bah, bph); //Compute outer angle. ba,bp
        auto gradba = rabn * (identity3() - tensor_product(bah, bah));
        auto grad_a = rabn * (bph - bah * cost);
        auto P1 = tensor_product(bph - bah * cost, -bah * rabn*rabn);    
        auto P2 = -rabn * (tensor_product(bah, grad_a) + cost * gradba);
        return d_get(c.f_outer_angle_p,j) * harmonic_energy_hessian(d_get(c.outer_angle_p0,j), cost, grad_a, grad_a, P1+P2); //Harmonic Energy Hessian
    }

    INLINE mat3 outer_angle_hessian_p_b(const Constants& c) const{
        auto cost = dot(bah, bph); //Compute outer angle. ba,bp
        auto gradba = rabn * (tensor_product(bah, bah) - identity3());
        auto gradbp = rbpn * (tensor_product(bph, bph) - identity3());
        auto grad_b = rabn * (bah * cost - bph) + rbpn * (bph * cost - bah);
        auto grad_a = rabn * (bph - bah * cost);
        auto P1 = tensor_product(bph - bah * cost, bah * rabn*rabn);
        auto P3 = rabn * (gradbp - (tensor_product(bah, grad_b) + cost * gradba));
        return d_get(c.f_outer_angle_p,j) * harmonic_energy_hessian(d_get(c.outer_angle_p0,j), cost, grad_a, grad_b, P1 + P3); //Harmonic Energy Hessian
    }

    INLINE mat3 outer_angle_hessian_p_p(const Constants& c) const{
        auto cost = dot(bah, bph); //Compute outer angle. ba,bp
        auto gradbp = rbpn * (identity3() - tensor_product(bph, bph));
        auto grad_a = rabn * (bph - bah * cost);
        auto grad_p = rbpn * (bah - bph * cost);
        auto P1 = rabn * (gradbp - tensor_product(bah, grad_p));
        return d_get(c.f_outer_angle_p,j) * harmonic_energy_hessian(d_get(c.outer_angle_p0,j), cost, grad_a, grad_p, P1); //Harmonic Energy Hessian
    }   

    INLINE auto dihedral_hessian_terms(const Constants& c) const{
        auto cbh = -bch;
        auto cost1 = dot(abh, cbh);
        auto cost2 = dot(cbh, dbh);
        auto sint1 = SQRT(1 - cost1*cost1);
        auto sint2 = SQRT(1 - cost2*cost2);
        auto cot1 = cost1/sint1;
        auto csc1 = device_real_t(1.)/sint1;
        auto csc2 = device_real_t(1.)/sint2;
        auto nabc = cross(abh, cbh) * csc1;
        auto nbcd = cross(dbh, cbh) * csc2;
        auto cosb = dot(nabc, nbcd);
        auto Coeff = cosb * csc1 * rabn;
        auto F1 = abh * sint1;
        auto F2 = cross(cbh, nbcd) / cosb;
        auto F3 = cot1 * (abh * cost1 - cbh);
        auto F = F1 - F2 + F3;
        auto GradACosb = cosb * rabn * csc1 * (abh * sint1 - cross(cbh, nbcd) / cosb + cot1 * (abh * cost1 - cbh)); 
        return std::tuple{cbh, cost1, cost2, sint1, sint2, cot1, csc1, csc2, nabc, nbcd, cosb, Coeff, F, GradACosb};
    }

    // $\nabla_a(\nabla_a(\cos(\theta)))$
    INLINE mat3 dihedral_hessian_a(const Constants& c) const{
        auto [cbh, cost1, cost2, sint1, sint2, cot1, csc1, csc2, nabc, nbcd, cosb, Coeff, F, GradACosb] = dihedral_hessian_terms(c);
        auto GradARab = abh * rabn * rabn;
        auto GradAabh = (tensor_product(abh, abh) - identity3()) * rabn;
        auto GradASint1 = -(abh * cost1 - cbh) * cost1 * rabn * csc1;
        auto GradAcsc1 = -GradASint1 * csc1 * csc1;
        auto GradAcost1 = (abh * cost1 - cbh) * rabn;
        auto GradAcot1 = (sint1 * GradAcost1 - cost1 * GradASint1) * csc1 * csc1;
        auto GradACoeff = GradACosb * rabn * csc1 + cosb * (GradARab * csc1 + GradAcsc1 * rabn);
        auto GradAF1 = GradAabh * sint1 + tensor_product(abh, GradASint1);
        auto GradAF2 = tensor_product(cross(cbh, nbcd), -GradACosb / (cosb * cosb));
        auto GradAF3 = tensor_product(abh * cost1 - cbh, GradAcot1) + cot1 * (tensor_product(abh, GradAcost1) + cost1 * GradAabh);
        auto GradAF = GradAF1 - GradAF2 + GradAF3;
        auto GradAGradCosb = tensor_product(F, GradACoeff) + Coeff * GradAF;
        return d_get(c.f_inner_dihedral,j) * harmonic_energy_hessian(d_get(c.inner_dih0,j), cosb, GradACosb, GradACosb, GradAGradCosb); //Harmonic Energy Hessian

    }

    // $\nabla_b(\nabla_a(\cos(\beta)))$
    INLINE mat3 dihedral_hessian_b(const Constants& c) const{
        auto [cbh, cost1, cost2, sint1, sint2, cot1, csc1, csc2, nabc, nbcd, cosb, Coeff, F, GradACosb] = dihedral_hessian_terms(c);
        auto grad_b_sint1 = -((cbh - abh*cost1)*rabn + (abh - cbh*cost1)*rbcn)*cost1 * csc1;
        auto grad_b_sint2 = -((cbh - dbh*cost2)*rdbn + (dbh - cbh*cost2)*rbcn)*cost2 * csc2;
        auto grad_b_ab_cross_cb_dot_nbcd = (rbcn * (cross(nbcd, abh) - dot(nbcd, cross(abh, cbh))*cbh) - rabn * (cross(nbcd, cbh) - dot(nbcd, cross(cbh, abh))*abh)); 
        auto grad_b_db_cross_cb_dot_nabc = (rbcn * (cross(nabc, dbh) - dot(nabc, cross(dbh, cbh))*cbh) - rdbn * (cross(nabc, cbh) - dot(nabc, cross(cbh, dbh))*dbh));
        auto P1 = (grad_b_ab_cross_cb_dot_nbcd*sint1 - (dot(nbcd, cross(abh, cbh)))*grad_b_sint1)*csc1*csc1;
        auto P2 = (grad_b_db_cross_cb_dot_nabc*sint2 - (dot(nabc, cross(dbh, cbh)))*grad_b_sint2)*csc2*csc2;
        auto grad_b = P1 + P2;
        auto GradBRab = -abh*rabn*rabn;
        auto GradBabh = (identity3() - tensor_product(abh, abh))*rabn;
        auto GradBcbh = (identity3() - tensor_product(cbh, cbh))*rbcn;
        auto GradBdbh = (identity3() - tensor_product(dbh, dbh))*rdbn;
        auto GradBnbcd = ((cross(dbh, GradBcbh) - cross(cbh, GradBdbh))* sint2 - tensor_product(cross(dbh,cbh), grad_b_sint2))*csc2*csc2;
        auto GradBcsc1 = -grad_b_sint1 * csc1*csc1;
        auto GradBcost1 = (cbh - abh*cost1)*rabn + (abh - cbh*cost1)*rbcn;
        auto GradBcot1 = (sint1 * GradBcost1 - cost1 * grad_b_sint1)*csc1*csc1;
        auto GradBCoeff = grad_b*rabn*csc1 + cosb*(GradBRab * csc1 + GradBcsc1*rabn);
        auto GradBF1 = GradBabh*sint1 + tensor_product(abh, grad_b_sint1);
        auto GradBF2 = tensor_product(cross(cbh,nbcd), -grad_b/(cosb*cosb)) + (cross(cbh, GradBnbcd) - cross(nbcd, GradBcbh))/cosb;
        auto GradBF3 = tensor_product(abh*cost1-cbh, GradBcot1) + cot1*(tensor_product(abh,GradBcost1) + cost1*GradBabh - GradBcbh);
        auto GradBF = GradBF1 - GradBF2 + GradBF3;
        auto GradBGradCosb = tensor_product(F, GradBCoeff) + Coeff*GradBF;
        return d_get(c.f_inner_dihedral,j) * harmonic_energy_hessian(d_get(c.inner_dih0,j), cosb, GradACosb, grad_b, GradBGradCosb); //Harmonic Energy Hessian
    }

    // $\nabla_c(\nabla_a(\cos(\theta)))$
    INLINE mat3 dihedral_hessian_c(const Constants& c) const{
        auto [cbh, cost1, cost2, sint1, sint2, cot1, csc1, csc2, nabc, nbcd, cosb, Coeff, F, GradACosb] = dihedral_hessian_terms(c);
        auto grad_c_sint1 = -(cbh*cost1 - abh)*cost1*csc1*rbcn;
        auto grad_c_sint2 = -(cbh*cost2 - dbh)*cost2*csc2*rbcn;
        auto grad_c_ab_cross_cb_dot_nabc = rbcn * (dot(nabc, cross(dbh, cbh))*cbh - cross(nabc, dbh));
        auto grad_c_db_cross_cb_dot_nbcd = rbcn * (dot(nbcd, cross(abh, cbh))*cbh - cross(nbcd, abh));
        auto P1 = (grad_c_ab_cross_cb_dot_nabc*sint2 - (dot(nabc, cross(dbh, cbh)))*grad_c_sint2)*csc2*csc2;
        auto P2 = (grad_c_db_cross_cb_dot_nbcd*sint1 - (dot(nbcd, cross(abh, cbh)))*grad_c_sint1)*csc1*csc1;
        auto grad_c = P1 + P2;

        auto GradCcbh   = (tensor_product(cbh, cbh) - identity3())*rbcn;
        auto GradCcsc1  = -grad_c_sint1 * csc1*csc1;
        auto GradCcost1 = (cbh*cost1 - abh)*rbcn;
        auto GradCcot1  = (sint1 * GradCcost1 - cost1 * grad_c_sint1)*csc1*csc1;
        auto GradCnbcd  = (cross(dbh, GradCcbh)*sint2 - tensor_product(cross(dbh, cbh), grad_c_sint2))*csc2*csc2;
        auto GradCCoeff = grad_c * rabn * csc1 + cosb*(GradCcsc1*rabn);
        auto GradCF1    = tensor_product(abh, grad_c_sint1);
        auto GradCF2    = tensor_product(cross(cbh, nbcd), -grad_c/(cosb*cosb)) + (cross(cbh, GradCnbcd) - cross(nbcd, GradCcbh))/cosb;
        auto GradCF3    = tensor_product(abh*cost1-cbh, GradCcot1) + cot1*(tensor_product(abh,GradCcost1) - GradCcbh);
        auto GradCF     = GradCF1 - GradCF2 + GradCF3;
        auto GradCGradCosb = tensor_product(F, GradCCoeff) + Coeff*GradCF;
        return d_get(c.f_inner_dihedral,j) * harmonic_energy_hessian(d_get(c.inner_dih0,j), cosb, GradACosb, grad_c, GradCGradCosb); //Harmonic Energy Hessian
    }

    // $\nabla_d(\nabla_a(\cos(\theta)))$
    INLINE mat3 dihedral_hessian_d(const Constants& c) const{
        auto [cbh, cost1, cost2, sint1, sint2, cot1, csc1, csc2, nabc, nbcd, cosb, Coeff, F, GradACosb] = dihedral_hessian_terms(c);
        auto GradDSint2 = -(dbh*cost2 - cbh)*cost2*csc2*rdbn;
        auto GradDDbCrossCbDotNabc = -rdbn * (dot(nabc, cross(cbh, dbh))*dbh - cross(nabc, cbh));
        auto grad_d = (GradDDbCrossCbDotNabc*sint2 - (dot(nabc, cross(dbh, cbh)))*GradDSint2)*csc2*csc2;

        auto GradDdbh = (tensor_product(dbh, dbh) - identity3())*rdbn;
        auto GradDnbcd = (-cross(cbh, GradDdbh)*sint2 - tensor_product(cross(dbh, cbh), GradDSint2))*csc2*csc2;
        auto GradDCoeff = grad_d * rabn * csc1;
        auto GradDF2 = tensor_product(cross(cbh, nbcd), -grad_d/(cosb*cosb)) + cross(cbh, GradDnbcd)/cosb;
        auto GradDF = - GradDF2;
        auto GradDGradCosb = tensor_product(F, GradDCoeff) + Coeff*GradDF;
        return d_get(c.f_inner_dihedral,j) * harmonic_energy_hessian(d_get(c.inner_dih0,j), cosb, GradACosb, grad_d, GradDGradCosb); //Harmonic Energy Hessian
    }

    INLINE auto outer_dihedral_hessian_a_terms(const Constants& c) const{
        auto mah = -amh;
        auto cosa = dot(abh, amh);
        auto cosm = dot(mah, mph);
        auto sina = sqrt(1 - cosa*cosa);
        auto sinm = sqrt(1 - cosm*cosm);
        auto cota = cosa/sina;
        auto cotm = cosm/sinm;
        auto csca = real_t(1.)/sina;
        auto cscm = real_t(1.)/sinm;
        auto nbam = cross(abh, amh)*csca;
        auto namp = cross(mah, mph)*cscm;
        auto cosb = dot(nbam, namp);
        auto F1 = abh*cosb;
        auto F2 = cross(amh, namp)*csca;
        auto F3 = amh*cosb;
        auto F4 = cross(namp, abh)*csca;
        auto G1 = abh*cosa * rabn;
        auto G2 = amh * rabn;
        auto G3 = amh*cosa * ramn;
        auto G4 = abh * ramn;
        auto H1 = cross(mph, nbam);
        auto H2 = mah*cosb*sinm;
        auto H3 = cotm*cosb*(mph - mah*cosm);
        auto C1 = cota*cosb*csca;
        auto C2 = ramn * cscm;
        auto GradAcosb = (F1 - F2)*rabn + (F3 - F4)*ramn + C1*(G1 - G2 + G3 - G4) + C2*(H1 - H2 + H3);
        return std::tuple(cosa, cosm, sina, sinm, cota, cotm, csca, cscm, nbam, namp, cosb, F1, F2, F3, F4, G1, G2, G3, G4, H1, H2, H3, C1, C2, GradAcosb);
    }

    INLINE mat3 outer_dihedral_hessian_a_a(const Constants& c) const{
        auto [cosa, cosm, sina, sinm, cota, cotm, csca, cscm, nbam, namp, cosb, F1, F2, F3, F4, G1, G2, G3, G4, H1, H2, H3, C1, C2, GradAcosb] = outer_dihedral_hessian_a_terms(c);
        auto mah = -amh;

        auto GradAcosa = (amh*cosa - abh)*ramn + (abh*cosa - amh)*rabn;
        auto GradAsina = -cosa*csca * GradAcosa;
        auto GradAcota = (sina * GradAcosa - cosa * GradAsina)*csca*csca;
        auto GradAcsca = -GradAsina*csca*csca;
        auto GradAcosm = (mph - mah*cosm)*ramn;
        auto GradAsinm = -cosm*cscm * GradAcosm;
        auto GradAcotm = (sinm * GradAcosm - cosm * GradAsinm)*cscm*cscm;
        auto GradAcscm = -GradAsinm*cscm*cscm;
        
        auto GradAab = (tensor_product(abh, abh) - identity3())*rabn;
        auto GradArabn = abh*rabn*rabn;
        auto GradAam = (tensor_product(amh, amh) - identity3())*ramn;
        auto GradAramn = amh*ramn*ramn;
        auto GradAma = (identity3() - tensor_product(mah,mah))*ramn;
        auto GradAnbam = ((cross(abh, GradAam) - cross(amh, GradAab))*sina  - tensor_product(cross(abh, amh), GradAsina))*csca*csca;
        auto GradAnamp = (( - cross(mph, GradAma))*sinm - tensor_product(cross(mah, mph), GradAsinm))*cscm*cscm;
        auto GradAF1 = tensor_product(abh, GradAcosb) + cosb*GradAab;
        auto GradAF2 = ((cross(amh, GradAnamp) - cross(namp, GradAam))*sina - tensor_product(cross(amh, namp), GradAsina))*csca*csca;
        auto GradAF3 = tensor_product(amh, GradAcosb) + cosb*GradAam;
        auto GradAF4 = ((cross(namp, GradAab) - cross(abh, GradAnamp))*sina - tensor_product(cross(namp, abh), GradAsina))*csca*csca;

        auto GradAG1 = tensor_product(abh, (GradAcosa*rabn + GradArabn*cosa)) + cosa*rabn * GradAab;
        auto GradAG2 = tensor_product(amh, GradArabn) + GradAam*rabn;
        auto GradAG3 = tensor_product(amh, (GradAcosa*ramn + GradAramn*cosa)) + cosa*ramn * GradAam;
        auto GradAG4 = tensor_product(abh, GradAramn) + GradAab*ramn;

        auto GradAH1 = cross(mph,GradAnbam);
        auto GradAH2 = tensor_product(mah, (cosb*GradAsinm + GradAcosb*sinm)) + cosb*sinm*GradAma;
        auto GradAH3 = tensor_product(mph - mah*cosm, (GradAcotm*cosb + GradAcosb*cotm)) + cotm*cosb*(- (GradAma*cosm + tensor_product(mah,GradAcosm)));

        auto GradAC1 = GradAcota * cosb*csca + cota* (GradAcosb*csca + cosb*GradAcsca);
        auto GradAC2 = GradAramn*cscm + GradAcscm*ramn;

        auto GradAF = rabn * (GradAF1 - GradAF2)  + tensor_product(F1 - F2,GradArabn) + ramn * (GradAF3 - GradAF4) + tensor_product(F3 - F4,GradAramn);
        auto GradAG = C1 * (GradAG1 - GradAG2 + GradAG3 - GradAG4) + tensor_product(G1 - G2 + G3 - G4, GradAC1);
        auto GradAH = C2 * (GradAH1 - GradAH2 + GradAH3) + tensor_product(H1 - H2 + H3, GradAC2);

        auto GradGradAcosb = GradAF + GradAG + GradAH;
        return d_get(c.f_outer_dihedral,j) * harmonic_energy_hessian(d_get(c.outer_dih0_a,j), cosb, GradAcosb, GradAcosb, GradGradAcosb); //Harmonic Energy Hessian
    }

    INLINE mat3 outer_dihedral_hessian_a_b(const Constants& c) const{
        auto [cosa, cosm, sina, sinm, cota, cotm, csca, cscm, nbam, namp, cosb, F1, F2, F3, F4, G1, G2, G3, G4, H1, H2, H3, C1, C2, GradAcosb] = outer_dihedral_hessian_a_terms(c);
        auto mah = -amh;

        auto GradBab = (identity3() - tensor_product(abh, abh))*rabn;
        auto GradBrabn = -abh*rabn*rabn;
        auto GradBcosa = (amh - cosa*abh)*rabn;
        auto GradBsina = -cosa/sina * GradBcosa;
        auto GradBcota = (sina * GradBcosa - cosa * GradBsina)*csca*csca;
        auto GradBcsca = -GradBsina*csca*csca;
        auto GradBnbam = (( - cross(amh, GradBab))*sina  - tensor_product(cross(abh, amh), GradBsina))*csca*csca;
        auto GradBcosb = (-(cross(namp,amh) - dot(namp,cross(amh,abh))*abh)*sina*rabn - dot(namp,cross(abh,amh))*GradBsina)*csca*csca;

        auto GradBF1 = tensor_product(abh, GradBcosb) + cosb*GradBab;
        auto GradBF2 = (- tensor_product(cross(amh, namp), GradBsina))*csca*csca;
        auto GradBF3 = tensor_product(amh, GradBcosb);
        auto GradBF4 = ((cross(namp, GradBab) )*sina - tensor_product(cross(namp, abh), GradBsina))*csca*csca;

        auto GradBG1 = tensor_product(abh, (GradBcosa*rabn + GradBrabn*cosa)) + cosa*rabn * GradBab;
        auto GradBG2 = tensor_product(amh, GradBrabn);
        auto GradBG3 = tensor_product(amh, (GradBcosa*ramn));
        auto GradBG4 = GradBab*ramn;

        auto GradBH1 = cross(mph, GradBnbam);
        auto GradBH2 = tensor_product(mah, GradBcosb*sinm);
        auto GradBH3 = tensor_product(mph - mah*cosm, (GradBcosb*cotm));

        auto GradBC1 = GradBcota * cosb*csca + cota* (GradBcosb * csca + cosb*GradBcsca);

        auto GradBF = rabn * (GradBF1 - GradBF2)  + tensor_product(F1 - F2,GradBrabn) + ramn * (GradBF3 - GradBF4);
        auto GradBG = C1 * (GradBG1 - GradBG2 + GradBG3 - GradBG4) + tensor_product(G1 - G2 + G3 - G4, GradBC1);
        auto GradBH = C2 * (GradBH1 - GradBH2 + GradBH3);

        auto GradGradBcosb = GradBF + GradBG + GradBH;
        return d_get(c.f_outer_dihedral,j) * harmonic_energy_hessian(d_get(c.outer_dih0_a,j), cosb, GradAcosb, GradBcosb, GradGradBcosb); //Harmonic Energy Hessian
    }

    INLINE mat3 outer_dihedral_hessian_a_m(const Constants& c) const{
        auto [cosa, cosm, sina, sinm, cota, cotm, csca, cscm, nbam, namp, cosb, F1, F2, F3, F4, G1, G2, G3, G4, H1, H2, H3, C1, C2, GradAcosb] = outer_dihedral_hessian_a_terms(c);
        auto mah = -amh;

        auto GradMcosa = (abh - cosa*amh)*ramn;
        auto GradMcosm = (mah*cosm - mph)*ramn + (mph*cosm - mah)*rmpn;
        auto GradMsina = -cosa/sina * GradMcosa;
        auto GradMsinm = -cosm/sinm * GradMcosm;

        auto GradMam = (identity3() - tensor_product(amh, amh)) * ramn;
        auto GradMramn = -amh * ramn * ramn;
        auto GradMma = (tensor_product(amh, amh) - identity3()) * ramn;
        auto GradMmp = (tensor_product(mph, mph) - identity3()) * rmpn;
        auto GradMcota = (sina * GradMcosa - cosa * GradMsina) * csca * csca;
        auto GradMcotm = (sinm * GradMcosm - cosm * GradMsinm) * cscm * cscm;
        auto GradMcsca = -GradMsina * csca * csca;
        auto GradMcscm = -GradMsinm * cscm * cscm;
        auto GradMnbam = ((cross(abh, GradMam))*sina - tensor_product(cross(abh, amh), GradMsina)) * csca * csca;
        auto GradMnamp = ((cross(mah, GradMmp) - cross(mph, GradMma))*sinm - tensor_product(cross(mah, mph), GradMsinm)) * cscm * cscm;
       
        auto cosbP1 = (((cross(namp,abh) - dot(namp,cross(abh,amh))*amh)*sina*ramn - dot(namp,cross(abh,amh))*GradMsina)*csca*csca);
        auto cosbP2 = (((rmpn*cscm)*(dot(nbam,cross(mah, mph))*mph - cross(nbam,mah)) - (cscm*ramn)*(dot(nbam,cross(mph, mah))*mah - cross(nbam,mph))) - dot(nbam,cross(mah, mph))*GradMsinm*cscm*cscm);
        auto GradMcosb = cosbP1 + cosbP2;   
        auto GradMF1 = tensor_product(abh, GradMcosb);
        auto GradMF2 = ((cross(amh, GradMnamp) - cross(namp, GradMam))*sina - tensor_product(cross(amh, namp), GradMsina)) * csca * csca;
        auto GradMF3 = tensor_product(amh, GradMcosb) + cosb*GradMam;
        auto GradMF4 = (( - cross(abh, GradMnamp))*sina - tensor_product(cross(namp, abh), GradMsina)) * csca * csca;

        auto GradMG1 = tensor_product(abh, GradMcosa*rabn);
        auto GradMG2 = GradMam*rabn;
        auto GradMG3 = tensor_product(amh, GradMcosa*ramn + GradMramn*cosa) + cosa*ramn*GradMam;
        auto GradMG4 = tensor_product(abh, GradMramn);

        auto GradMH1 = cross(mph, GradMnbam) - cross(nbam, GradMmp);
        auto GradMH2 = tensor_product(mah, (cosb*GradMsinm + GradMcosb*sinm)) + cosb*sinm*GradMma;
        auto GradMH3 = tensor_product(mph - mah*cosm, (GradMcotm*cosb + GradMcosb*cotm)) + cotm*cosb*(GradMmp - (GradMma*cosm + tensor_product(mah,GradMcosm)));

        auto GradMC1 = GradMcota * cosb * csca + cota* (GradMcosb*csca + cosb*GradMcsca);
        auto GradMC2 = GradMramn * cscm + GradMcscm*ramn;

        auto GradMF = rabn * (GradMF1 - GradMF2) + ramn * (GradMF3 - GradMF4) + tensor_product(F3 - F4, GradMramn);
        auto GradMG = C1 * (GradMG1 - GradMG2 + GradMG3 - GradMG4) + tensor_product(G1 - G2 + G3 - G4, GradMC1);
        auto GradMH = C2 * (GradMH1 - GradMH2 + GradMH3) + tensor_product(H1 - H2 + H3, GradMC2);

        auto GradGradMcosb = GradMF + GradMG + GradMH;
        return d_get(c.f_outer_dihedral,j) * harmonic_energy_hessian(d_get(c.outer_dih0_a,j), cosb, GradAcosb, GradMcosb, GradGradMcosb); //Harmonic Energy Hessian
    }

    INLINE mat3 outer_dihedral_hessian_a_p(const Constants& c) const{
        auto [cosa, cosm, sina, sinm, cota, cotm, csca, cscm, nbam, namp, cosb, F1, F2, F3, F4, G1, G2, G3, G4, H1, H2, H3, C1, C2, GradAcosb] = outer_dihedral_hessian_a_terms(c);
        auto mah = -amh;

        auto GradPcosm = (mah - mph * cosm) * rmpn;
        auto GradPsinm = -cosm * cscm * GradPcosm;

        auto GradPmp = (identity3() - tensor_product(mph, mph)) * rmpn;
        auto GradPcotm = (sinm * GradPcosm - cosm * GradPsinm) * cscm*cscm;
        auto GradPcscm = -GradPsinm * cscm*cscm;
        auto GradPnbam = -tensor_product(cross(abh, amh), GradAcosb) * csca*csca;
        auto GradPnamp = (cross(mah, GradPmp) * sinm - tensor_product(cross(mah, mph), GradPsinm)) * cscm*cscm;
        auto GradPcosb = (sinm * rmpn * (cross(nbam, mah) - dot(nbam,cross(mah, mph))*mph) - dot(nbam,cross(mah, mph)) * GradPsinm) * cscm*cscm;

        auto GradPF1 = tensor_product(abh, GradPcosb);
        auto GradPF2 = ((cross(amh,GradPnamp))*sina ) * csca*csca;
        auto GradPF3 = tensor_product(amh, GradPcosb);
        auto GradPF4 = ((- cross(abh,GradPnamp))*sina ) * csca*csca;

        auto GradPH1 =  - cross(nbam, GradPmp);
        auto GradPH2 = tensor_product(mah, (cosb*GradPsinm + GradPcosb*sinm));
        auto GradPH3 = tensor_product(mph - mah*cosm, (GradPcotm*cosb + GradPcosb*cotm)) + cotm*cosb*(GradPmp - (tensor_product(mah,GradPcosm)));

        auto GradPC1 = cota* (GradPcosb * csca);
        auto GradPC2 = GradPcscm * ramn;

        auto GradPF = (GradPF1 - GradPF2) * rabn + (GradPF3 - GradPF4) * ramn;
        auto GradPG = tensor_product(G1 - G2 + G3 - G4, GradPC1);
        auto GradPH = C2 * (GradPH1 - GradPH2 + GradPH3) + tensor_product(H1 - H2 + H3, GradPC2);

        auto GradGradPcosb = GradPF + GradPG + GradPH;
        return d_get(c.f_outer_dihedral,j) * harmonic_energy_hessian(d_get(c.outer_dih0_a,j), cosb, GradAcosb, GradPcosb, GradGradPcosb); //Harmonic Energy Hessian
    }

    INLINE auto outer_dihedral_hessian_m_terms() const {
        auto pmh = -mph;
        auto cosm = dot(mbh, mph);
        auto cosp = dot(pmh, pah);
        auto sinm = sqrt(1 - cosm*cosm);
        auto sinp = sqrt(1 - cosp*cosp);
        auto cscm = 1/sinm;
        auto cscp = 1/sinp;
        auto cotp = cosp/sinp;
        auto cotm = cosm/sinm;
        auto nbmp = cross(mbh, mph)*cscm;
        auto nmpa = cross(pmh, pah)*cscp;
        auto cosb = dot(nbmp, nmpa);
        auto F = cross(nbmp, pmh) * rapn * cscp;
        auto G = pah*cosb * rapn;
        auto K1 = cotp*cosb;
        auto K2 = rapn*cscp;
        auto K = K1 * K2;
        auto H = K * (pmh - pah*cosp);
        auto GradAcosb = F - G + H;
        return std::tuple(cosm, cosp, sinm, sinp, cscm, cscp, cotp, cotm, nbmp, nmpa, cosb, F, G, H, K1, K2, K, GradAcosb);
    }

    INLINE mat3 outer_dihedral_hessian_m_a(const Constants& c) const{
        auto [cosm, cosp, sinm, sinp, cscm, cscp, cotp, cotm, nbmp, nmpa, cosb, F, G, H, K1, K2, K, GradAcosb] = outer_dihedral_hessian_m_terms();
        auto pmh = -mph;

        auto GradAcosp = (pmh - pah*cosp)*rapn;
        auto GradAsinp = -cosp*cscp * GradAcosp;
        auto GradAcscp = -GradAsinp*cscp*cscp;
        auto GradAcotp = (sinp * GradAcosp - cosp * GradAsinp)*cscp*cscp;

        auto GradApah = (identity3() - tensor_product(pah, pah))*rapn;
        auto GradArpan = -pah*rapn*rapn;
        auto GradAF = tensor_product(cross(nbmp,pmh), GradArpan * cscp + GradAcscp * rapn);
        auto GradAG = tensor_product(pah, GradAcosb*rapn + GradArpan*cosb) + GradApah * cosb * rapn;
        auto GradAK1 = GradAcotp * cosb + cotp * GradAcosb;
        auto GradAK2 = GradArpan * cscp + GradAcscp * rapn;
        auto GradAK = GradAK1 * K2 + K1 * GradAK2;
        auto GradAH = tensor_product(pmh-pah*cosp, GradAK) + K * (- GradApah * cosp - tensor_product(pah, GradAcosp));

        auto GradGradAcosb = GradAF - GradAG + GradAH;
        return d_get(c.f_outer_dihedral,j) * harmonic_energy_hessian(d_get(c.outer_dih0_m,j), cosb, GradAcosb, GradAcosb, GradGradAcosb); //Harmonic Energy Hessian
    }

    INLINE mat3 outer_dihedral_hessian_m_b(const Constants& c) const{
        auto [cosm, cosp, sinm, sinp, cscm, cscp, cotp, cotm, nbmp, nmpa, cosb, F, G, H, K1, K2, K, GradAcosb] = outer_dihedral_hessian_m_terms();
        auto pmh = -mph;

        auto GradBcosm = (mph - mbh*cosm)*rbmn;
        auto GradBsinm = -cosm/sinm * GradBcosm;
        auto GradBcscm = -GradBsinm*cscm*cscm;
        
        auto GradBcosb = - (cross(nmpa, mph) - dot(cross(nmpa, mph), mbh)*mbh)*rbmn*cscm - dot(nmpa, cross(mbh, mph)) * GradBsinm*cscm*cscm;
        
        auto GradBmbh = (identity3() - tensor_product(mbh, mbh))*rbmn;
        auto GradBnbmp = (-cross(mph, GradBmbh)*sinm - tensor_product(cross(mbh, mph), GradBsinm))*cscm*cscm;
        auto GradBF = -cross(pmh, GradBnbmp) * rapn * cscp;
        auto GradBG = tensor_product(pah, GradBcosb*rapn);
        auto GradBK1 = cotp * GradBcosb;
        auto GradBK = GradBK1 * K2;
        auto GradBH = tensor_product(pmh - pah*cosp, GradBK);
        auto GradGradBcosb = GradBF - GradBG + GradBH;
        return d_get(c.f_outer_dihedral,j) * harmonic_energy_hessian(d_get(c.outer_dih0_m,j), cosb, GradAcosb, GradBcosb, GradGradBcosb); //Harmonic Energy Hessian
    }

    INLINE mat3 outer_dihedral_hessian_m_m(const Constants& c) const {
        auto [cosm, cosp, sinm, sinp, cscm, cscp, cotp, cotm, nbmp, nmpa, cosb, F, G, H, K1, K2, K, GradAcosb] = outer_dihedral_hessian_m_terms();
        auto pmh = -mph;

        auto GradMcosm = (mbh*cosm - mph)*rbmn + (mph*cosm - mbh)*rmpn;
        auto GradMsinm = -cosm/sinm * GradMcosm;
        auto GradMcscm = -GradMsinm*cscm*cscm;
        auto GradMcosp = (pah - pmh*cosp)*rmpn;
        auto GradMsinp = -cosp * cscp * GradMcosp;
        auto GradMcscp = -GradMsinp*cscp*cscp;
        auto GradMcotp = (sinp * GradMcosp - cosp * GradMsinp)*cscp*cscp;

        auto GradMmph = (tensor_product(mph,mph) - identity3())*rmpn;
        auto GradMpmh = (identity3() - tensor_product(pmh,pmh))*rmpn;
        auto GradMmbh = (tensor_product(mbh,mbh) - identity3())*rbmn;

        auto GradMnbmp = ((cross(mbh, GradMmph) - cross(mph, GradMmbh))*sinm - tensor_product(cross(mbh,mph),GradMsinm))*cscm*cscm;
        auto GradMcosbP1 = (((dot(cross(nmpa, mbh),mph)*mph - cross(nmpa, mbh))*rmpn - (dot(cross(nmpa, mph),mbh)*mbh - cross(nmpa, mph))*rbmn)*sinm - dot(nmpa,cross(mbh,mph))*GradMsinm)*cscm*cscm;
        auto GradMcosbP2 = ((- (cross(nbmp, pah) - dot(cross(nbmp, pah),pmh)*pmh)*rmpn)*sinp - dot(nbmp,cross(pmh,pah))*GradMsinp)*cscp*cscp;
        auto GradMcosb = GradMcosbP1 + GradMcosbP2;

        auto GradMF = (cross(nbmp, GradMpmh) - cross(pmh, GradMnbmp))*rapn*cscp + tensor_product(cross(nbmp,pmh),GradMcscp*rapn);
        auto GradMG = tensor_product(pah, GradMcosb*rapn);
        auto GradMK1 = GradMcotp * cosb + cotp * GradMcosb;
        auto GradMK2 = GradMcscp * rapn;
        auto GradMK = GradMK1 * K2 + K1 * GradMK2;
        auto GradMH = tensor_product((pmh - pah*cosp), GradMK) + K * (GradMpmh - tensor_product(pah,GradMcosp));
        auto GradGradMcosb = GradMF - GradMG + GradMH;
        return d_get(c.f_outer_dihedral,j) * harmonic_energy_hessian(d_get(c.outer_dih0_m,j), cosb, GradAcosb, GradMcosb, GradGradMcosb); //Harmonic Energy Hessian
    }

    INLINE mat3 outer_dihedral_hessian_m_p(const Constants& c) const {
        auto [cosm, cosp, sinm, sinp, cscm, cscp, cotp, cotm, nbmp, nmpa, cosb, F, G, H, K1, K2, K, GradAcosb] = outer_dihedral_hessian_m_terms();
        auto pmh = -mph;

        auto GradPcosm = (mbh - mph*cosm)*rmpn;
        auto GradPsinm = -cosm/sinm * GradPcosm;
        auto GradPcscm = -GradPsinm * cscm*cscm;
        auto GradPcosp = (pmh*cosp - pah)*rmpn + (pah*cosp - pmh)*rapn;
        auto GradPsinp = -cosp*cscp * GradPcosp;
        auto GradPcscp = -GradPsinp * cscp*cscp;
        auto GradPcotp = (sinp * GradPcosp - cosp * GradPsinp)*cscp*cscp;

        auto GradPmph = (identity3() - tensor_product(mph,mph))*rmpn;
        auto GradPpmh = (tensor_product(pmh,pmh) - identity3())*rmpn;
        auto GradPpah = (tensor_product(pah,pah) - identity3())*rapn;

        auto GradPrpan = pah*rapn*rapn;

        auto GradPnbmp = ((cross(mbh,GradPmph))*sinm - tensor_product(cross(mbh,mph),GradPsinm))*cscm*cscm;
        auto GradPcosbP1 = (((cross(nmpa, mbh) - dot(cross(nmpa, mbh),mph)*mph)*rmpn)*sinm - dot(nmpa,cross(mbh,mph))*GradPsinm)*cscm*cscm;
        auto GradPcosbP2 = ( ((dot(cross(nbmp, pmh),pah)*pah - cross(nbmp, pmh))*rapn  - (dot(cross(nbmp, pah),pmh)*pmh - cross(nbmp, pah))*rmpn) * sinp - dot(nbmp,cross(pmh,pah))*GradPsinp)* cscp*cscp;
        auto GradPcosb = GradPcosbP1 + GradPcosbP2;

        auto GradPF = (cross(nbmp, GradPpmh) - cross(pmh, GradPnbmp))*rapn*cscp + tensor_product(cross(nbmp,pmh),GradPcscp*rapn + GradPrpan*cscp);
        auto GradPG = tensor_product(pah, GradPcosb*rapn + GradPrpan*cosb) +  GradPpah * rapn*cosb;
        auto GradPK1 = GradPcotp * cosb + cotp * GradPcosb;
        auto GradPK2 = GradPrpan * cscp + GradPcscp * rapn;
        auto GradPK = GradPK1 * K2 + K1 * GradPK2;
        auto GradPH = tensor_product((pmh - pah*cosp), GradPK) + K * (GradPpmh - tensor_product(pah,GradPcosp) - GradPpah * cosp);

        auto GradGradPcosb = GradPF - GradPG + GradPH;
        return d_get(c.f_outer_dihedral,j) * harmonic_energy_hessian(d_get(c.outer_dih0_m,j), cosb, GradAcosb, GradPcosb, GradGradPcosb); //Harmonic Energy Hessian
    }
    
    INLINE auto outer_dihedral_hessian_p_terms() const{
        auto pah = -aph;
        auto cosa = dot(aph,amh);
        auto cosp = dot(pbh,pah);
        auto sina = sqrt(1 - cosa*cosa);
        auto sinp = sqrt(1 - cosp*cosp);
        auto csca = 1/sina;
        auto cscp = 1/sinp;
        auto cotp = cosp/sinp;
        auto cota = cosa/sina;
        auto nbpa = cross(pbh,pah)*cscp;
        auto npam = cross(aph,amh)*csca;
        auto cosb = dot(nbpa,npam);

        auto rpan = rapn;
        auto C1 = rpan*cscp;
        auto C2 = cota*cosb*csca;
        auto F1 = cross(npam,pbh);
        auto F2 = pah*cosb*sinp;
        auto F3 = cotp*cosb*(pbh - pah*cosp);
        auto G1 = aph*cosb;
        auto G2 = cross(amh,nbpa)*csca;
        auto G3 = amh*cosb;
        auto G4 = cross(nbpa,aph)*csca;
        auto H1 = aph*cosa * rapn;
        auto H2 = amh * rapn;
        auto H3 = amh*cosa * ramn;
        auto H4 = aph * ramn;

        auto GradAcosb = C1 * (F1- F2 + F3) + rapn * (G1 - G2)  + ramn * (G3 - G4) + C2 * (H1 - H2 + H3 - H4);
        return std::tuple(cosa, cosp, sina, sinp, csca, cscp, cotp, cota, nbpa, npam, cosb, C1, C2, F1, F2, F3, G1, G2, G3, G4, H1, H2, H3, H4, GradAcosb);
    }

    INLINE mat3 outer_dihedral_hessian_p_a(const Constants& c) const {
        auto [cosa, cosp, sina, sinp, csca, cscp, cotp, cota, nbpa, npam, cosb, C1, C2, F1, F2, F3, G1, G2, G3, G4, H1, H2, H3, H4, GradAcosb] = outer_dihedral_hessian_p_terms();

        auto GradAcosa = (aph*cosa - amh)*rapn + (amh*cosa - aph)*ramn;
        auto GradAsina = -cosa * csca * GradAcosa;
        auto GradAcsca = -GradAsina * csca * csca;
        auto GradAcota = (sina * GradAcosa - cosa * GradAsina)*csca*csca;
        auto GradAcosp = (pbh - pah*cosp)*rapn;
        auto GradAsinp = -cosp * cscp * GradAcosp;
        auto GradAcscp = -GradAsinp * cscp * cscp;
        auto GradAcotp = (sinp * GradAcosp - cosp * GradAsinp)*cscp*cscp;

        auto GradAamh = (tensor_product(amh,amh) - identity3())*ramn;
        auto GradAramn = amh*ramn*ramn;
        auto GradAaph = (tensor_product(aph,aph) - identity3())*rapn;
        auto GradArapn = aph*rapn*rapn;

        auto GradApah = (identity3() - tensor_product(pah,pah))*rapn;
        auto GradArpan = -pah*rapn*rapn;

        auto GradAC1 = cscp * GradArpan + rapn * GradAcscp;
        auto GradAC2 = cosb * (cota* GradAcsca + csca * GradAcota) + cota * csca * GradAcosb;

        auto GradAnpam = (sina*(cross(aph, GradAamh) - cross(amh, GradAaph)) - tensor_product(cross(aph,amh),GradAsina))*csca*csca;

        auto GradAF1 = - cross(pbh,GradAnpam);
        auto GradAF2 = tensor_product(pah, sinp*GradAcosb + cosb*GradAsinp) + GradApah *sinp*cosb;
        auto GradAF3 = cotp*cosb*( - tensor_product(pah,GradAcosp) - GradApah*cosp) + tensor_product(pbh - pah*cosp,GradAcotp*cosb) + tensor_product(pbh - pah*cosp,cotp*GradAcosb);

        auto GradAnbpa = (sinp*(cross(pbh, GradApah)) - tensor_product(cross(pbh,pah),GradAsinp))*cscp*cscp;

        auto GradAG1 = tensor_product(aph,GradAcosb) + GradAaph*cosb;
        auto GradAG2 = tensor_product(cross(amh,nbpa), GradAcsca) + csca*(cross(amh,GradAnbpa) - cross(nbpa,GradAamh));
        auto GradAG3 = tensor_product(amh,GradAcosb) + GradAamh*cosb;
        auto GradAG4 = tensor_product(cross(nbpa,aph), GradAcsca) + csca*(cross(nbpa,GradAaph) - cross(aph,GradAnbpa));

        auto GradAH1 = tensor_product(aph,GradAcosa*rapn + GradArapn*cosa) + GradAaph*cosa*rapn;
        auto GradAH2 = tensor_product(amh,GradArapn) + GradAamh*rapn;
        auto GradAH3 = tensor_product(amh,GradAcosa*ramn + GradAramn*cosa) + GradAamh*cosa*ramn;
        auto GradAH4 = tensor_product(aph,GradAramn) + GradAaph*ramn;

        auto GradGradAcosb = C1 * (GradAF1 - GradAF2 + GradAF3) + tensor_product(F1 - F2 + F3, GradAC1) + rapn * (GradAG1 - GradAG2) + tensor_product(G1 - G2,GradArapn) + ramn * (GradAG3 - GradAG4) + tensor_product(G3 - G4,GradAramn) + C2 * (GradAH1 - GradAH2 + GradAH3 - GradAH4) + tensor_product(H1 - H2 + H3 - H4,GradAC2);
        return d_get(c.f_outer_dihedral,j) * harmonic_energy_hessian(d_get(c.outer_dih0_p,j), cosb, GradAcosb, GradAcosb, GradGradAcosb); //Harmonic Energy Hessian
    }

    INLINE mat3 outer_dihedral_hessian_p_b(const Constants& c) const {
        auto [cosa, cosp, sina, sinp, csca, cscp, cotp, cota, nbpa, npam, cosb, C1, C2, F1, F2, F3, G1, G2, G3, G4, H1, H2, H3, H4, GradAcosb] = outer_dihedral_hessian_p_terms();

        auto GradBcosp = (pah - pbh*cosp)*rbpn;
        auto GradBsinp = -cosp * cscp * GradBcosp;
        auto GradBcscp = -GradBsinp * cscp * cscp;
        auto GradBcotp = (sinp * GradBcosp - cosp * GradBsinp)*cscp*cscp;

        auto GradBpbh = (identity3() - tensor_product(pbh,pbh))*rbpn;
        auto GradBrpbn = -pbh*rbpn*rbpn;
        auto GradBcosb = ((dot(npam,cross(pah,pbh))*pbh - cross(npam,pah))*sinp*rbpn - dot(npam,cross(pbh,pah))*GradBsinp)*cscp*cscp;

        auto GradBC1 = GradBcscp * rapn;
        auto GradBC2 = cota * csca * GradBcosb;

        auto GradBF1 = cross(npam,GradBpbh);
        auto GradBF2 = tensor_product(pah, sinp*GradBcosb  + cosb*GradBsinp);
        auto GradBF3 = cotp*cosb*(GradBpbh - tensor_product(pah,GradBcosp)) + tensor_product(pbh - pah*cosp,GradBcotp*cosb) + tensor_product(pbh - pah*cosp,cotp*GradBcosb);

        auto GradBnbpa = (sinp*( - cross(pah, GradBpbh)) - tensor_product(cross(pbh,pah),GradBsinp))*cscp*cscp;
        auto GradBG1 = tensor_product(aph,GradBcosb);
        auto GradBG2 =  csca*(cross(amh,GradBnbpa));
        auto GradBG3 = tensor_product(amh,GradBcosb);
        auto GradBG4 = csca*(- cross(aph,GradBnbpa));

        auto GradGradBcosb = C1 * (GradBF1 - GradBF2 + GradBF3) + tensor_product(F1 - F2 + F3, GradBC1) + rapn * (GradBG1 - GradBG2) + ramn * (GradBG3 - GradBG4) + tensor_product(H1 - H2 + H3 - H4,GradBC2);

        return d_get(c.f_outer_dihedral,j) * harmonic_energy_hessian(d_get(c.outer_dih0_p,j), cosb, GradAcosb, GradBcosb, GradGradBcosb); //Harmonic Energy Hessian
    }

    INLINE mat3 outer_dihedral_hessian_p_m(const Constants& c) const{
        auto [cosa, cosp, sina, sinp, csca, cscp, cotp, cota, nbpa, npam, cosb, C1, C2, F1, F2, F3, G1, G2, G3, G4, H1, H2, H3, H4, GradAcosb] = outer_dihedral_hessian_p_terms();

        auto GradMcosa = (aph - amh*cosa)*ramn;
        auto GradMsina = -cosa * csca * GradMcosa;
        auto GradMcsca = -GradMsina * csca * csca;
        auto GradMcota = (sina * GradMcosa - cosa * GradMsina)*csca*csca;
        auto GradMcosb = ((cross(nbpa,aph) - dot(cross(nbpa,aph),amh)*amh)*sina*ramn - dot(nbpa,cross(aph,amh))*GradMsina)*csca*csca;

        auto GradMamh = (identity3() - tensor_product(amh,amh))*ramn;
        auto GradMramn = -amh*ramn*ramn;

        auto GradMC2 = cosb * (cota* GradMcsca + csca * GradMcota) + cota * csca * GradMcosb;
        auto GradMnpam = (sina*(cross(aph, GradMamh)) - tensor_product(cross(aph,amh),GradMsina))*csca*csca;
        auto GradMF1 = - cross(pbh,GradMnpam);
        auto GradMF2 = tensor_product(pah, sinp*GradMcosb);
        auto GradMF3 = tensor_product(pbh - pah*cosp,cotp*GradMcosb);

        auto GradMG1 = tensor_product(aph,GradMcosb);
        auto GradMG2 = tensor_product(cross(amh,nbpa), GradMcsca) + csca*(- cross(nbpa,GradMamh));
        auto GradMG3 = tensor_product(amh,GradMcosb) + GradMamh*cosb;
        auto GradMG4 = tensor_product(cross(nbpa,aph), GradMcsca);

        auto GradMH1 = tensor_product(aph,GradMcosa*rapn);
        auto GradMH2 = GradMamh*rapn;
        auto GradMH3 = tensor_product(amh,GradMcosa*ramn + GradMramn*cosa) + GradMamh*cosa*ramn;
        auto GradMH4 = tensor_product(aph,GradMramn);

        auto GradGradMcosb = C1 * (GradMF1 - GradMF2 + GradMF3) + rapn * (GradMG1 - GradMG2) + ramn * (GradMG3 - GradMG4) + tensor_product(G3 - G4,GradMramn) + C2 * (GradMH1 - GradMH2 + GradMH3 - GradMH4) + tensor_product(H1 - H2 + H3 - H4, GradMC2);
        return d_get(c.f_outer_dihedral,j) * harmonic_energy_hessian(d_get(c.outer_dih0_p,j), cosb, GradAcosb, GradMcosb, GradGradMcosb); //Harmonic Energy Hessian
    }

    INLINE mat3 outer_dihedral_hessian_p_p(const Constants& c) const{
        auto [cosa, cosp, sina, sinp, csca, cscp, cotp, cota, nbpa, npam, cosb, C1, C2, F1, F2, F3, G1, G2, G3, G4, H1, H2, H3, H4, GradAcosb] = outer_dihedral_hessian_p_terms();

        auto GradPcosa = (amh - aph*cosa)*rapn;
        auto GradPsina = -cosa * csca * GradPcosa;
        auto GradPcsca = -GradPsina * csca * csca;
        auto GradPcota = (sina * GradPcosa - cosa * GradPsina)*csca*csca;
        auto GradPcosp = (pbh*cosp - pah)*rbpn + (pah*cosp - pbh)*rapn;
        auto GradPsinp = -cosp * cscp * GradPcosp;
        auto GradPcscp = -GradPsinp * cscp * cscp;
        auto GradPcotp = (sinp * GradPcosp - cosp * GradPsinp)*cscp*cscp;

        auto GradPaph = (identity3() - tensor_product(aph,aph))*rapn;
        auto GradPrapn = -aph*rapn*rapn;
        auto GradPpbh = (tensor_product(pbh,pbh) - identity3())*rbpn;
        auto GradPrpbn = pbh*rbpn*rbpn;
        auto GradPpah = (tensor_product(pah,pah) - identity3())*rapn;
        auto GradPrpan = pah*rapn*rapn;
        auto GradPcosbP1 = (((dot(cross(npam,pbh),pah)*pah - cross(npam,pbh))*rapn - (dot(cross(npam,pah),pbh)*pbh - cross(npam,pah))*rbpn)*sinp - dot(npam,cross(pbh,pah))*GradPsinp)*cscp*cscp;
        auto GradPcosbP2 = (-(cross(nbpa,amh) - dot(cross(nbpa,amh),aph)*aph)*sina*rapn - dot(nbpa,cross(aph,amh))*GradPsina)*csca*csca;
        auto GradPcosb = GradPcosbP1 + GradPcosbP2;

        auto GradPC1 = cscp * GradPrpan + rapn * GradPcscp;
        auto GradPC2 = cosb * (cota* GradPcsca + csca * GradPcota) + cota * csca * GradPcosb;

        auto GradPnpam = (sina*( -cross(amh,GradPaph)) - tensor_product(cross(aph,amh),GradPsina))*csca*csca;
        auto GradPF1 = cross(npam,GradPpbh) - cross(pbh,GradPnpam);
        auto GradPF2 = tensor_product(pah, sinp*GradPcosb + cosb*GradPsinp) + GradPpah *sinp*cosb;
        auto GradPF3 = cotp*cosb*(GradPpbh - tensor_product(pah,GradPcosp) - GradPpah*cosp) + tensor_product(pbh - pah*cosp,GradPcotp*cosb) + tensor_product(pbh - pah*cosp,cotp*GradPcosb);

        auto GradPG1 = tensor_product(aph,GradPcosb) + GradPaph*cosb;
        auto GradPnbpa = (sinp*(cross(pbh,GradPpah) - cross(pah,GradPpbh)) - tensor_product(cross(pbh,pah),GradPsinp))*cscp*cscp;
        auto GradPG2 = tensor_product(cross(amh,nbpa), GradPcsca) + csca*(cross(amh,GradPnbpa) );
        auto GradPG3 = tensor_product(amh,GradPcosb);
        auto GradPG4 = tensor_product(cross(nbpa,aph), GradPcsca) + csca*(cross(nbpa,GradPaph) - cross(aph,GradPnbpa));

        auto GradPH1 = tensor_product(aph,GradPcosa*rapn + GradPrapn*cosa) + GradPaph*cosa*rapn;
        auto GradPH2 = tensor_product(amh,GradPrapn) ;
        auto GradPH3 = tensor_product(amh,GradPcosa*ramn);
        auto GradPH4 = GradPaph*ramn;
    
        auto GradGradPcosb = C1 * (GradPF1 - GradPF2 + GradPF3) + tensor_product(F1 - F2 + F3, GradPC1) + rapn * (GradPG1 - GradPG2) + tensor_product(G1 - G2,GradPrapn) + ramn * (GradPG3 - GradPG4) + C2 * (GradPH1 - GradPH2 + GradPH3 - GradPH4) + tensor_product(H1 - H2 + H3 - H4, GradPC2);
        return d_get(c.f_outer_dihedral,j) * harmonic_energy_hessian(d_get(c.outer_dih0_p,j), cosb, GradAcosb, GradPcosb, GradGradPcosb); //Harmonic Energy Hessian
    }

    //Computes gradient related to bending of outer angles. ~20 FLOPs
    /**
     * @brief Compute the gradient of the outer angle-m term.
     * @param c The constants for the threadIdx^th node.
     * @return The gradient of the outer angle-m term.
    */
    INLINE coord3d outer_angle_gradient_m(const Constants& c) const
    {
        real_t cos_angle = -dot(abh, bmh); //Compute outer angle. ab,bm
        coord3d grad = (bmh + abh * cos_angle) * rabn; //Derivative of outer angles Eq. 30. Buster Thesis
        return d_get(c.f_outer_angle_m,j) * harmonic_energy_gradient(d_get(c.outer_angle_m0,j),cos_angle,grad); //Harmonic Energy Gradient: Eq. 30 multiplied by harmonic term.
    }

    /**
     * @brief Compute the gradient of the outer angle-p term.
     * @param c The constants for the threadIdx^th node.
     * @return The gradient of the outer angle-p term.
    */
    INLINE coord3d outer_angle_gradient_p(const Constants& c) const
    {   
        real_t cos_angle = -dot(abh, bph); //Compute outer angle. ab,bp
        coord3d grad = (bph + abh * cos_angle) * rabn; //Derivative of outer angles Eq. 28. Buster Thesis
        return d_get(c.f_outer_angle_p,j) * harmonic_energy_gradient(d_get(c.outer_angle_p0,j),cos_angle,grad); //Harmonic Energy Gradient: Eq. 28 multiplied by harmonic term.
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
        cos_b = dot(bah,bch); r_sin_b = (real_t)1.0/SQRT((real_t)1.0 - cos_b*cos_b); nabc = cross(bah, bch) * r_sin_b;
        cos_c = dot(-bch,cdh); r_sin_c = (real_t)1.0/SQRT((real_t)1.0 - cos_c*cos_c); nbcd = cross(-bch,cdh) * r_sin_c;

        real_t cos_beta = dot(nabc, nbcd); //Inner dihedral angle from planes abc,bcd.
        real_t cot_b = cos_b * r_sin_b * r_sin_b; //cos(b)/sin(b)^2

        //Derivative w.r.t. inner dihedral angle F and G in Eq. 26
        coord3d grad = cross(bch, nbcd) * r_sin_b * rabn - bah * cos_beta * rabn + (cot_b * cos_beta * rabn) * (bch - bah * cos_b);
        return d_get(c.f_inner_dihedral,j) * harmonic_energy_gradient(d_get(c.inner_dih0,j), cos_beta, grad); //Eq. 26.
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

        cos_a = dot(abh,amh); r_sin_a = (real_t)1.0/SQRT((real_t)1.0 - cos_a*cos_a); nbam_hat = cross(abh,amh) * r_sin_a;
        cos_m = dot(-amh,mph); r_sin_m = (real_t)1.0/SQRT((real_t)1.0 - cos_m*cos_m); namp_hat = cross(-amh,mph) * r_sin_m;
        
        real_t cos_beta = dot(nbam_hat, namp_hat); //Outer Dihedral angle bam, amp
        real_t cot_a = cos_a * r_sin_a * r_sin_a;
        real_t cot_m = cos_m * r_sin_m * r_sin_m;

        //Derivative w.r.t. outer dihedral angle, factorized version of Eq. 31.
        coord3d grad = cross(mph,nbam_hat)*ramn*r_sin_m - (cross(namp_hat,abh)*ramn + cross(amh,namp_hat)*rabn)*r_sin_a +
                        cos_beta*(abh*rabn + ramn * ((real_t)2.0*amh + cot_m*(mph+cos_m*amh)) - cot_a*(ramn*(abh - amh*cos_a) + rabn*(amh-abh*cos_a)));
        
        //Eq. 31 multiplied by harmonic term.
        return d_get(c.f_outer_dihedral,j) * harmonic_energy_gradient(d_get(c.outer_dih0_a,j), cos_beta, grad);
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
        cos_m = dot(mbh,mph);  r_sin_m = (real_t)1.0/SQRT((real_t)1.0 - cos_m*cos_m); nbmp_hat = cross(mbh,mph) * r_sin_m;
        cos_p = dot(-mph,pah); r_sin_p = (real_t)1.0/SQRT((real_t)1.0 - cos_p*cos_p); nmpa_hat = cross(-mph,pah) * r_sin_p;
        
        //Cosine to the outer dihedral angle constituted by the planes bmp and mpa
        real_t cos_beta = dot(nbmp_hat, nmpa_hat); //Outer dihedral angle bmp,mpa.
        real_t cot_p = cos_p * r_sin_p * r_sin_p;
        
        //Derivative w.r.t. outer dihedral angle, factorized version of Eq. 32.
        coord3d grad = rapn * (cot_p*cos_beta * (-mph - pah*cos_p) - cross(nbmp_hat, mph)*r_sin_p - pah*cos_beta );

        //Eq. 32 multiplied by harmonic term.
        return d_get(c.f_outer_dihedral,j) * harmonic_energy_gradient(d_get(c.outer_dih0_m,j), cos_beta, grad);
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
        cos_a = dot(aph,amh);  r_sin_a = (real_t)1.0/SQRT((real_t)1.0 - cos_a*cos_a); npam_hat = cross(aph,amh)  * r_sin_a;
        cos_p = dot(pbh,-aph); r_sin_p = (real_t)1.0/SQRT((real_t)1.0 - cos_p*cos_p); nbpa_hat = cross(pbh,-aph) * r_sin_p;

        real_t cos_beta = dot(nbpa_hat, npam_hat); //Outer dihedral angle bpa, pam.
        real_t cot_p = cos_p * r_sin_p * r_sin_p;
        real_t cot_a = cos_a * r_sin_a * r_sin_a;

        //Derivative w.r.t. outer dihedral angle, factorized version of Eq. 33.
        coord3d grad = cross(npam_hat,pbh)*rapn*r_sin_p - (cross(amh,nbpa_hat)*rapn + cross(nbpa_hat,aph)*ramn)*r_sin_a +
                        cos_beta*(amh*ramn + rapn * ((real_t)2.0*aph + cot_p*(pbh+cos_p*aph)) - cot_a*(rapn*(amh - aph*cos_a) + ramn*(aph-amh*cos_a)));
        
        //Eq. 33 multiplied by harmonic term.
        return d_get(c.f_outer_dihedral,j) * harmonic_energy_gradient(d_get(c.outer_dih0_p,j), cos_beta, grad);
    }

    // Internal coordinate gradients
    /**
     * @brief Compute the gradient of the bond length term.
     * @param c The constants for the threadIdx^th node.
     * @return The gradient of the bond length term.
    */
    INLINE coord3d bond_length_gradient(const Constants& c) const {
        return d_get(c.f_bond,j) * harmonic_energy_gradient(bond(),d_get(c.r0,j),abh); 
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
        case PEDERSEN:
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
        return (real_t)0.5 *d_get(c.f_bond,j) *harmonic_energy(bond(),d_get(c.r0,j));
    }
    /**
     * @brief Compute the total energy contribution of the bending terms.
     * @param c The constants for the threadIdx^th node.
     * @return The energy contribution of the bending terms.
    */
    INLINE real_t bend_energy(const Constants& c) const {
        return d_get(c.f_inner_angle,j)* harmonic_energy(angle(),d_get(c.angle0,j));
    }

    /**
     * @brief Compute the total energy contribution of the dihedral terms.
     * @param c The constants for the threadIdx^th node.
     * @return The energy contribution of the dihedral terms.
    */
    INLINE real_t dihedral_energy(const Constants& c) const {
        return d_get(c.f_inner_dihedral,j)* harmonic_energy(dihedral(),d_get(c.inner_dih0,j));
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
            break;
        case BOND:
            return bond_length_gradient(c);
            break;
        case ANGLE:
            return inner_angle_gradient(c);
            break;
        case DIH:
            return inner_dihedral_gradient(c);
            break;
        case ANGLE_M:
            return outer_angle_gradient_m(c);
            break;
        case ANGLE_P:
            return outer_angle_gradient_p(c);
            break;
        case DIH_A:
            return outer_dihedral_gradient_a(c);
            break;
        case DIH_M:
            return outer_dihedral_gradient_m(c);
            break;
        case DIH_P:
            return outer_dihedral_gradient_p(c);
            break;
        default:
            return bond_length_gradient(c) + angle_gradient(c) + dihedral_gradient(c);
            break;
        }
        return bond_length_gradient(c) + angle_gradient(c) + dihedral_gradient(c);
    }

    INLINE mat3 hessian_a(const Constants& c) const {
        switch (T)
        {
        case BOND:
            return bond_hessian_a(c);
            break;
        case ANGLE:
            return inner_angle_hessian_a(c);
            break;
        case ANGLE_M:
            return outer_angle_hessian_m_a(c);
            break;
        case ANGLE_P:
            return outer_angle_hessian_p_a(c);
            break;
        case DIH:
            return dihedral_hessian_a(c);
            break;
        case DIH_A:
            return outer_dihedral_hessian_a_a(c);
            break;
        case DIH_M:
            return outer_dihedral_hessian_m_a(c);
            break;
        case DIH_P:
            return outer_dihedral_hessian_p_a(c);
            break;
        default:
            return bond_hessian_a(c) + inner_angle_hessian_a(c) + outer_angle_hessian_m_a(c) + outer_angle_hessian_p_a(c) + dihedral_hessian_a(c) + outer_dihedral_hessian_a_a(c) + outer_dihedral_hessian_m_a(c) + outer_dihedral_hessian_p_a(c);
            break;
        }
    }

    INLINE mat3 hessian_b(const Constants& c) const {
        switch (T)
        {
        case BOND:
            return bond_hessian_b(c);
            break;
        case ANGLE:
            return inner_angle_hessian_b(c);
            break;
        case ANGLE_M:
            return outer_angle_hessian_m_b(c);
            break;
        case ANGLE_P:
            return outer_angle_hessian_p_b(c);
            break;
        case DIH:
            return dihedral_hessian_b(c);
            break;
        case DIH_A:
            return outer_dihedral_hessian_a_b(c);
            break;
        case DIH_M:
            return outer_dihedral_hessian_m_b(c);
            break;
        case DIH_P:
            return outer_dihedral_hessian_p_b(c);
            break;
        default:
            return bond_hessian_b(c) + inner_angle_hessian_b(c) + outer_angle_hessian_m_b(c) + outer_angle_hessian_p_b(c) + dihedral_hessian_b(c) + outer_dihedral_hessian_a_b(c) + outer_dihedral_hessian_m_b(c) + outer_dihedral_hessian_p_b(c);
            break;
        }
    }

    INLINE mat3 hessian_c(const Constants& c) const {
        switch (T){
        case BOND:
            return mat3();
            break;
        case ANGLE:
            return inner_angle_hessian_c(c);
            break;
        case ANGLE_M:
            return mat3();
            break;
        case ANGLE_P:
            return mat3();
            break;
        case DIH:
            return dihedral_hessian_c(c);
            break;
        case DIH_A:
            return mat3();
            break;
        case DIH_M:
            return mat3();
            break;
        case DIH_P:
            return mat3();
            break;
        default:
            return inner_angle_hessian_c(c) + dihedral_hessian_c(c);
            break;
        }
        return inner_angle_hessian_c(c) + dihedral_hessian_c(c);
    }

    INLINE mat3 hessian_d(const Constants& c) const{
        switch (T){
        case BOND:
            return mat3();
            break;
        case ANGLE:
            return mat3();
            break;
        case ANGLE_M:
            return mat3();
            break;
        case ANGLE_P:
            return mat3();
            break;
        case DIH:
            return dihedral_hessian_d(c);
            break;
        case DIH_A:
            return mat3();
            break;
        case DIH_M:
            return mat3();
            break;
        case DIH_P:
            return mat3();
            break;
        default:
            return dihedral_hessian_d(c);
            break;
        }
        return dihedral_hessian_d(c);
    }

    INLINE mat3 hessian_m(const Constants& c) const {
        switch(T)
        {
        case BOND:
            return mat3();
            break;
        case ANGLE:
            return mat3();
            break;
        case ANGLE_M:
            return outer_angle_hessian_m_m(c);
            break;
        case ANGLE_P:
            return mat3();
            break;
        case DIH:
            return mat3();
            break;
        case DIH_A:
            return outer_dihedral_hessian_a_m(c);
            break;
        case DIH_M:
            return outer_dihedral_hessian_m_m(c);
            break;
        case DIH_P:
            return outer_dihedral_hessian_p_m(c);
            break;
        default:
            return outer_angle_hessian_m_m(c) + outer_dihedral_hessian_a_m(c) + outer_dihedral_hessian_m_m(c) + outer_dihedral_hessian_p_m(c);
            break;
        }
        return outer_angle_hessian_m_m(c) + outer_dihedral_hessian_a_m(c) + outer_dihedral_hessian_m_m(c) + outer_dihedral_hessian_p_m(c);
    }

    INLINE mat3 hessian_p(const Constants& c) const {
        switch (T)
        {
        case BOND:
            return mat3();
            break;
        case ANGLE:
            return mat3();
            break;
        case ANGLE_M:
            return mat3();
            break;
        case ANGLE_P:
            return outer_angle_hessian_p_p(c);
            break;
        case DIH:
            return mat3();
            break;
        case DIH_A:
            return outer_dihedral_hessian_a_p(c);
            break;
        case DIH_M:
            return outer_dihedral_hessian_m_p(c);
            break;
        case DIH_P:
            return outer_dihedral_hessian_p_p(c);
            break;
        default:
            return outer_angle_hessian_p_p(c) + outer_dihedral_hessian_a_p(c) + outer_dihedral_hessian_m_p(c) + outer_dihedral_hessian_p_p(c);
            break;
        }
        return outer_angle_hessian_p_p(c) + outer_dihedral_hessian_a_p(c) + outer_dihedral_hessian_m_p(c) + outer_dihedral_hessian_p_p(c);
    }


    //Reciprocal lengths of arcs ab, ac, am, ap.
    real_t
        rab,
        rabn,
        racn,
        radn,
        ramn,
        rbmn,
        rbpn,
        rbcn,
        rdbn,
        rmpn,
        rapn;

    //Base Arcs,
    coord3d
        ab,
        ac,
        ad;

    //All normalized arcs required to perform energy & gradient calculations.
    //Note that all these arcs are cyclical the arc ab becomes: ab->ac->ad,  the arc ac becomes: ac->ad->ab , the arc bc becomes: bc->cd->db (For iterations 0, 1, 2)
    //As such the naming convention here is related to the arcs as they are used in the 0th iteration.
    coord3d 
        abh,
        ach,
        adh,
        bph,
        bmh,
        amh,
        aph,
        bah,
        bch,
        cdh,
        dbh,
        mph,
        mbh,
        pah,
        pbh;

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
INLINE coord3d gradient(const coord3d* X) const {
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

INLINE hessian_t hessian(coord3d* X) const {
    BLOCK_SYNC
    hessian_t hess(node_graph);
    for (uint8_t j = 0; j < 3; j++ ){
        ArcData arc = ArcData(j, X, node_graph);
        hess.A[0] += arc.hessian_a(constants);
        hess.A[1 + j] += arc.hessian_b(constants);
        hess.A[1 + (j + 1) % 3] += arc.hessian_c(constants);
        hess.A[1 + (j + 2) % 3] += arc.hessian_d(constants);
        hess.A[4 + j] += arc.hessian_m(constants);
        hess.A[7 + j] += arc.hessian_p(constants);
    }
    return hess;
}

//Uses finite difference to compute the hessian
INLINE hessian_t fd_hessian(coord3d* X, const float reldelta = 1e-7) const{
    hessian_t hess_fd(node_graph);
    for (uint16_t i = 0; i < blockDim.x; i++){
        for (uint8_t j = 0; j < 10; j++){
            auto node = hess_fd.indices[j];
            coord3d X0 = X[node];
            for (uint8_t k = 0; k < 3; k++){
                if (i == threadIdx.x){ X[node][k] = X0[k] + X0[k]*reldelta;}
                coord3d grad_X0_p = gradient(X);
                BLOCK_SYNC
                if (i == threadIdx.x){ X[node][k] = X0[k] - X0[k]*reldelta;} 
                coord3d grad_X0_m = gradient(X);
                BLOCK_SYNC
                if (i == threadIdx.x){ 
                    hess_fd.A[j][0][k] = (grad_X0_p[0] - grad_X0_m[0])/(2*X0[k]*reldelta);
                    hess_fd.A[j][1][k] = (grad_X0_p[1] - grad_X0_m[1])/(2*X0[k]*reldelta);
                    hess_fd.A[j][2][k] = (grad_X0_p[2] - grad_X0_m[2])/(2*X0[k]*reldelta);
                    X[node][k] = X0[k];
                }
                BLOCK_SYNC
            }
        }
    }
    return hess_fd;
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
};

template <ForcefieldType T> __global__ void compute_hessians_(const IsomerBatch B, CuArray<device_real_t> Hess, CuArray<device_node_t> Cols){
    DEVICE_TYPEDEFS;
    extern __shared__ real_t smem[];
    clear_cache(smem,Block_Size_Pow_2);
    auto limit = ((B.isomer_capacity + gridDim.x - 1) / gridDim.x ) * gridDim.x;  //Fast ceiling integer division.
    for (int isomer_idx = blockIdx.x; isomer_idx < limit; isomer_idx += gridDim.x){
    BLOCK_SYNC
    if (isomer_idx < B.isomer_capacity)
      if(B.statuses[isomer_idx] == IsomerStatus::CONVERGED) { //Avoid illegal memory access
    size_t offset = isomer_idx * blockDim.x;
    Constants constants     = Constants(B, isomer_idx);
    NodeNeighbours node_graph    = NodeNeighbours(B, isomer_idx, smem);
    ForceField<T> FF           = ForceField<T>(node_graph, constants, smem);
    coord3d* X              = reinterpret_cast<coord3d*>(smem) + B.n_atoms;
    assign(X[threadIdx.x],reinterpret_cast<std::array<float,3>*>(B.X+offset*3)[threadIdx.x]);
    BLOCK_SYNC
    auto hess = FF.hessian(X);
    int n_cols = 10*3;
    int n_rows = blockDim.x*3;
    int hess_stride = n_cols*n_rows;
    int toff = isomer_idx*hess_stride + threadIdx.x*n_cols*3;
    for (size_t i = 0; i < 3; i++) //rows
    for (size_t j = 0; j < 10; j++) //cols / 3
    {   
        for (size_t k = 0; k < 3; k++)
        {
            Cols.data[toff + i*n_cols + j*3 + k] = hess.indices[j]*3 + k;
            Hess.data[toff + i*n_cols + j*3 + k] = hess.A[j][i][k];
        }    
    }
    //auto result = hess.lanczos_iteration<3>(reinterpret_cast<real_t*>(smem));
    //auto lambda = hess.power_iteration(reinterpret_cast<coord3d*>(smem));
    
    }}
}

template <ForcefieldType T> __global__ void compute_hessians_fd_(const IsomerBatch B, CuArray<device_real_t> Hess, CuArray<device_node_t> Cols, float reldelta){
    DEVICE_TYPEDEFS;
    extern __shared__ real_t smem[];
    clear_cache(smem,Block_Size_Pow_2);
    auto limit = ((B.isomer_capacity + gridDim.x - 1) / gridDim.x ) * gridDim.x;  //Fast ceiling integer division.
    for (int isomer_idx = blockIdx.x; isomer_idx < limit; isomer_idx += gridDim.x){
    BLOCK_SYNC
    if (isomer_idx < B.isomer_capacity) //Avoid illegal memory access
      if(B.statuses[isomer_idx] == IsomerStatus::CONVERGED) { // Only compute Hessian for valid geometries
    size_t offset = isomer_idx * blockDim.x;
    Constants constants     = Constants(B, isomer_idx);
    NodeNeighbours node_graph    = NodeNeighbours(B, isomer_idx, smem);
    ForceField<T> FF           = ForceField<T>(node_graph, constants, smem);
    coord3d* X              = reinterpret_cast<coord3d*>(smem + B.n_atoms);
    assign(X[threadIdx.x],reinterpret_cast<std::array<float,3>*>(B.X+offset*3)[threadIdx.x]);
    BLOCK_SYNC
    auto hess = FF.fd_hessian(X, reldelta);
    int n_cols = 10*3;
    int n_rows = blockDim.x*3;
    int hess_stride = n_cols*n_rows;
    int toff = isomer_idx*hess_stride + threadIdx.x*n_cols*3;
    for (size_t i = 0; i < 3; i++) //rows
    for (size_t j = 0; j < 10; j++) //cols / 3
    {
        for (size_t k = 0; k < 3; k++)
        {
            Cols.data[toff + i*n_cols + j*3 + k] = hess.indices[j]*3 + k;
            Hess.data[toff + i*n_cols + j*3 + k] = hess.A[j][i][k];
        }    
    }
    }}
}

float kernel_time = 0.0;

template <ForcefieldType T>
cudaError_t compute_hessians(const IsomerBatch& B, CuArray<device_real_t>& hessians, CuArray<device_node_t>& cols, const LaunchCtx& ctx, const LaunchPolicy policy){
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
    static LaunchDims dims((void*)compute_hessians_<T>, B.n_atoms, smem, B.isomer_capacity);
    dims.update_dims((void*)compute_hessians_<T>, B.n_atoms, smem, B.isomer_capacity);
    void* kargs[]{(void*)&B, (void*)&hessians, (void*)&cols};

    cudaEventRecord(start[dev], ctx.stream);
    auto error = safeCudaKernelCall((void*)compute_hessians_<T>, dims.get_grid(), dims.get_block(), kargs, smem, ctx.stream);
    cudaEventRecord(stop[dev], ctx.stream);
    
    if(policy == LaunchPolicy::SYNC) {
        ctx.wait();
        cudaEventElapsedTime(&single_kernel_time, start[dev], stop[dev]);
        kernel_time += single_kernel_time;
    }
    printLastCudaError("Hessian Computation Failed: ");
    first_call[dev] = false;
    return error;
}

template <ForcefieldType T>
cudaError_t compute_hessians_fd(const IsomerBatch& B, CuArray<device_real_t>& hessians, CuArray<device_node_t>& cols, const float reldelta, const LaunchCtx& ctx, const LaunchPolicy policy){
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
    static LaunchDims dims((void*)compute_hessians_fd_<T>, B.n_atoms, smem, B.isomer_capacity);
    dims.update_dims((void*)compute_hessians_fd_<T>, B.n_atoms, smem, B.isomer_capacity);
    void* kargs[]{(void*)&B, (void*)&hessians, (void*)&cols , (void*)&reldelta};

    cudaEventRecord(start[dev], ctx.stream);
    auto error = safeCudaKernelCall((void*)compute_hessians_fd_<T>, dims.get_grid(), dims.get_block(), kargs, smem, ctx.stream);
    cudaEventRecord(stop[dev], ctx.stream);
    
    if(policy == LaunchPolicy::SYNC) {
        ctx.wait();
        cudaEventElapsedTime(&single_kernel_time, start[dev], stop[dev]);
        kernel_time += single_kernel_time;
    }
    printLastCudaError("Finite Difference Hessian Computation Failed: ");
    first_call[dev] = false;
    return error;
}


void declaration(){
    IsomerBatch B(20,1,DEVICE_BUFFER);
    CuArray<float> arr(1);
    CuArray<device_real_t> hessians(1);
    CuArray<device_node_t> cols(1);

    compute_hessians<PEDERSEN>(B, hessians, cols);
    compute_hessians<BOND>(B, hessians, cols);
    compute_hessians<ANGLE>(B, hessians, cols);
    compute_hessians<ANGLE_M>(B, hessians, cols);
    compute_hessians<ANGLE_P>(B, hessians, cols);
    compute_hessians<DIH>(B, hessians, cols);
    compute_hessians<DIH_A>(B, hessians, cols);
    compute_hessians<DIH_M>(B, hessians, cols);
    compute_hessians<DIH_P>(B, hessians, cols);

    compute_hessians_fd<PEDERSEN>(B, hessians, cols, 0.0001);
    compute_hessians_fd<BOND>(B, hessians, cols,0.0001);
    compute_hessians_fd<ANGLE>(B, hessians, cols,0.0001);
    compute_hessians_fd<ANGLE_M>(B, hessians, cols,0.0001);
    compute_hessians_fd<ANGLE_P>(B, hessians, cols,0.0001);
    compute_hessians_fd<DIH>(B, hessians, cols,0.0001);
    compute_hessians_fd<DIH_A>(B, hessians, cols,0.0001);
    compute_hessians_fd<DIH_M>(B, hessians, cols,0.0001);
    compute_hessians_fd<DIH_P>(B, hessians, cols,0.0001);

}

} // namespace isomer
} // namespace timemachine
