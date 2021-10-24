

#include "cuda_runtime.h"
#include <cuda_runtime_api.h>
#include <cuda.h>
#include "device_launch_parameters.h"
#include <stdio.h>
#define getLastCudaError(x) 
#include <iostream>
#include <fstream>
#include <chrono>
#include <cuda_fp16.h>
#include <cuda_bf16.h>
#include "fullerenes/gpu/isomerspace_forcefield.hh"

#define BLOCK_SYNC cg::sync(cg::this_thread_block())
#define GRID_SYNC cg::sync(cg::this_grid())

typedef IsomerspaceForcefield::device_real_t device_real_t;
typedef IsomerspaceForcefield::device_node_t device_node_t;
typedef GPU_REAL3 device_coord3d;
typedef GPU_NODE3 device_node3;


#include "coord3d.cu"
#include "helper_functions.cu"
#include "data_structs.cu"

using namespace std::literals;
namespace cg = cooperative_groups;

/** This struct was made to reduce signature cluttering of device functions, it is simply a container for default arguments which are shared between functions**/
struct ForceField{
    typedef device_coord3d coord3d;
    typedef device_node3 node3;
    typedef device_node_t node_t;
    typedef device_real_t real_t;


    coord3d r1,                         //d_i+1 In conjugated gradient algorithm  Buster p. 36
            r0;                         //d_0 In conjugated gradient algorithm Buster p. 36

    const NodeGraph node_graph;         //Contains face-information and neighbour-information. Both of which are constant in the lifespan of this struct. 
    const Constants constants; //Contains force-constants and equillibrium-parameters. Constant in the lifespan of this struct.

    size_t node_id = threadIdx.x;
    real_t* sdata;                      //Pointer to start of L1 cache array, used exclusively for reduction.

    __device__ ForceField(  const NodeGraph &G,
                            const Constants &c, 
                            real_t* sdata): node_graph(G), constants(c), sdata(sdata) {}


//Container for all energy and gradient evaluations with respect to an arc, eg. AB, AC or AD.
struct ArcData{
    //124 FLOPs;

    __device__ ArcData(const u_char j, const coord3d* __restrict__ X, const NodeGraph& G){   
        this->j = j;   
        node_t a = threadIdx.x;
        real_t r_rmp;
        coord3d ap, am, ab, ac, ad, mp;
        coord3d X_a = X[a]; coord3d X_b = X[d_get(G.neighbours,j)];
        //printf("Index: %d \n", a*3 + j);

        //Compute the arcs ab, ac, ad, bp, bm, ap, am, mp, bc and cd
        ab = (X_b - X_a);  r_rab = bond_length(ab); ab_hat = r_rab * ab;
        ac = (X[d_get(G.neighbours,(j+1)%3)] - X_a); r_rac = bond_length(ac); ac_hat = r_rac * ac; rab = non_resciprocal_bond_length(ab);
        ad = (X[d_get(G.neighbours,(j+2)%3)] - X_a); r_rad = bond_length(ad); ad_hat = r_rad * ad;
        
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
    INLINE real_t harmonic_energy(const real_t p0, const real_t p) const{
        return (real_t)0.5*(p-p0)*(p-p0);
    }
    //4 FLOPs
    INLINE coord3d  harmonic_energy_gradient(const real_t p0, const real_t p, const coord3d gradp) const{
        return (p-p0)*gradp;     
    }
    //1 FLOP
    INLINE real_t bond() const {return (real_t)rab;}
    //5 FLOPs
    INLINE real_t angle() const {return dot(ab_hat,ac_hat);}

    //Returns outer angle m, used only diagnostically.
    INLINE real_t outer_angle_m() const {return -dot(ab_hat, bm_hat);} //Compute outer angle. ab,bm

    //Returns outer angle p, used only diagnostically.
    INLINE real_t outer_angle_p() const{return -dot(ab_hat, bp_hat);} //Compute outer angle. ab,bp

    //Returns the inner dihedral angle for the current arc. Used here only for energy calculation, 
    //otherwise embedded in dihedral computation because the planes and angles that make up the dihedral angle computation are required for derivative computation.
    //50 FLOPs
    INLINE real_t dihedral() const 
    { 
        coord3d nabc, nbcd; real_t cos_b, cos_c, r_sin_b, r_sin_c;
        cos_b = dot(ba_hat,bc_hat);  r_sin_b = (real_t)1.0/sqrt((real_t)1.0 - cos_b*cos_b); nabc = cross(ba_hat, bc_hat) * r_sin_b;
        cos_c = dot(-bc_hat,cd_hat); r_sin_c = (real_t)1.0/sqrt((real_t)1.0 - cos_c*cos_c); nbcd = cross(-bc_hat,cd_hat) * r_sin_c;
        return dot(nabc, nbcd);
    }
    //Returns the Outer-dihedral-a wrt. current arc, only accessed diagnostically (internal coordinate).
    INLINE real_t outer_dihedral_a() const
    {
        coord3d nbam_hat, namp_hat; real_t cos_a, cos_m, r_sin_a, r_sin_m;
        cos_a = dot(ab_hat,am_hat); r_sin_a = (real_t)1.0/sqrt((real_t)1.0 - cos_a*cos_a); nbam_hat = cross(ab_hat,am_hat) * r_sin_a;
        cos_m = dot(-am_hat,mp_hat); r_sin_m = (real_t)1.0/sqrt((real_t)1.0 - cos_m*cos_m); namp_hat = cross(-am_hat,mp_hat) * r_sin_m;
        real_t cos_beta = dot(nbam_hat, namp_hat); //Outer Dihedral angle bam, amp
        return cos_beta;
    }
    //Returns the Outer-dihedral-m wrt. current arc, only accessed diagnostically (internal coordinate).
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
    INLINE coord3d inner_angle_gradient(const Constants& c) const
    {   
        real_t cos_angle = angle(); //Inner angle of arcs ab,ac.
        coord3d grad = cos_angle * (ab_hat * r_rab + ac_hat * r_rac) - ab_hat * r_rac - ac_hat* r_rab; //Derivative of inner angle: Eq. 21. 
        return d_get(c.f_inner_angle,j) * harmonic_energy_gradient(d_get(c.angle0,j), cos_angle, grad); //Harmonic Energy Gradient: Eq. 21. multiplied by harmonic term.
    }
    //Computes gradient related to bending of outer angles. ~20 FLOPs
    INLINE coord3d outer_angle_gradient_m(const Constants& c) const
    {
        real_t cos_angle = -dot(ab_hat, bm_hat); //Compute outer angle. ab,bm
        coord3d grad = (bm_hat + ab_hat * cos_angle) * r_rab; //Derivative of outer angles Eq. 30. Buster Thesis
        return d_get(c.f_outer_angle_m,j) * harmonic_energy_gradient(d_get(c.outer_angle_m0,j),cos_angle,grad); //Harmonic Energy Gradient: Eq. 30 multiplied by harmonic term.
    }
    INLINE coord3d outer_angle_gradient_p(const Constants& c) const
    {
        real_t cos_angle = -dot(ab_hat, bp_hat); //Compute outer angle. ab,bp
        coord3d grad = (bp_hat + ab_hat * cos_angle) * r_rab; //Derivative of outer angles Eq. 28. Buster Thesis
        return d_get(c.f_outer_angle_p,j) * harmonic_energy_gradient(d_get(c.outer_angle_p0,j),cos_angle,grad); //Harmonic Energy Gradient: Eq. 28 multiplied by harmonic term.
    }
    // Chain rule terms for dihedral calculation
    //Computes gradient related to dihedral/out-of-plane term. ~75 FLOPs
    INLINE coord3d inner_dihedral_gradient(const Constants& c) const
    {
        coord3d nabc, nbcd; real_t cos_b, cos_c, r_sin_b, r_sin_c;
        cos_b = dot(ba_hat,bc_hat); r_sin_b = (real_t)1.0/sqrt((real_t)1.0 - cos_b*cos_b); nabc = cross(ba_hat, bc_hat) * r_sin_b;
        cos_c = dot(-bc_hat,cd_hat); r_sin_c = (real_t)1.0/sqrt((real_t)1.0 - cos_c*cos_c); nbcd = cross(-bc_hat,cd_hat) * r_sin_c;

        real_t cos_beta = dot(nabc, nbcd); //Inner dihedral angle from planes abc,bcd.
        real_t cot_b = cos_b * r_sin_b * r_sin_b; //cos(b)/sin(b)^2

        //Derivative w.r.t. inner dihedral angle F and G in Eq. 26
        coord3d grad = cross(bc_hat, nbcd) * r_sin_b * r_rab - ba_hat * cos_beta * r_rab + (cot_b * cos_beta * r_rab) * (bc_hat - ba_hat * cos_b);

        return d_get(c.f_inner_dihedral,j) * harmonic_energy_gradient(d_get(c.inner_dih0,j), cos_beta, grad); //Eq. 26.
    }

    //Computes gradient from dihedral angles constituted by the planes bam, amp ~162 FLOPs
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

        return d_get(c.f_outer_dihedral,j) * harmonic_energy_gradient(d_get(c.outer_dih0_a,j), cos_beta, grad);
    }

    //Computes gradient from dihedral angles constituted by the planes nbmp, nmpa ~92 FLOPs
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
        return d_get(c.f_outer_dihedral,j) * harmonic_energy_gradient(d_get(c.outer_dih0_m,j), cos_beta, grad);
    }

    //Computes gradient from dihedral angles constituted by the planes bpa, pam ~162 FLOPs
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
        return d_get(c.f_outer_dihedral,j) * harmonic_energy_gradient(d_get(c.outer_dih0_p,j), cos_beta, grad);
    }
    // Internal coordinate gradients
    INLINE coord3d bond_length_gradient(const Constants& c) const { return d_get(c.f_bond,j) * harmonic_energy_gradient(bond(),d_get(c.r0,j),ab_hat);}
    //Sum of angular gradient components.
    INLINE coord3d angle_gradient(const Constants& c) const { return inner_angle_gradient(c) + outer_angle_gradient_p(c) + outer_angle_gradient_m(c);}
    //Sum of inner and outer dihedral gradient components.
    INLINE coord3d dihedral_gradient(const Constants& c) const { return inner_dihedral_gradient(c) + outer_dihedral_gradient_a(c) + outer_dihedral_gradient_m(c) + outer_dihedral_gradient_p(c);}
    //coord3d flatness()             const { return ;  }   
    
    INLINE real_t bond_energy(const Constants& c) const {return (real_t)0.5 *d_get(c.f_bond,j) *harmonic_energy(bond(),d_get(c.r0,j));}
    INLINE real_t bend_energy(const Constants& c) const {return d_get(c.f_inner_angle,j)* harmonic_energy(angle(),d_get(c.angle0,j));}
    INLINE real_t dihedral_energy(const Constants& c) const {return d_get(c.f_inner_dihedral,j)* harmonic_energy(dihedral(),d_get(c.inner_dih0,j));}
    //Harmonic energy contribution from bond stretching, angular bending and dihedral angle bending.
    //71 FLOPs
    INLINE real_t energy(const Constants& c) const {return bond_energy(c) + bend_energy(c) + dihedral_energy(c); }
    //Sum of bond, angular and dihedral gradient components.
    INLINE coord3d gradient(const Constants& c) const{return bond_length_gradient(c) + angle_gradient(c) + dihedral_gradient(c);}

    
    uint8_t j;

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

    /*
    All normalized arcs required to perform energy & gradient calculations.
    Note that all these arcs are cyclical the arc ab becomes: ab->ac->ad,  the arc ac becomes: ac->ad->ab , the arc bc becomes: bc->cd->db (For iterations 0, 1, 2)
    As such the naming convention here is related to the arcs as they are used in the 0th iteration. */
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
};


INLINE coord3d gradient(coord3d* X) const {
    BLOCK_SYNC;
    coord3d grad = {0.0, 0.0, 0.0};

    for (uint8_t j = 0; j < 3; j++ ){
        ArcData arc = ArcData(j, X, node_graph);
        grad += arc.gradient(constants);
    }

    return grad;
}

INLINE real_t energy(coord3d* X) const {
    BLOCK_SYNC;
    real_t arc_energy = (real_t)0.0;

    //(71 + 124) * 3 * N  = 585*N FLOPs
    for (uint8_t j = 0; j < 3; j++ ){
        ArcData arc = ArcData(j, X, node_graph);
        arc_energy += arc.energy(constants);
    }
    return reduction(sdata, arc_energy);;

}

INLINE real_t gradnorm(coord3d* X, coord3d& d)const {
    return reduction(sdata, dot(-gradient(X),d));
}

//Bracketing method designed to find upper bound for linesearch method that matches 
//reference python implementation by Buster.
INLINE real_t FindLineSearchBound(coord3d* X, coord3d& r0, coord3d* X1){
    real_t bound = 1e-5;
    bool negative_grad = true;
    size_t iter = 0;
    while (negative_grad && iter < 30)
    {   
        bound *= (real_t)1.5;
        X1[node_id] = X[node_id] + bound * r0;
        real_t gradsum = reduction(sdata, dot(gradient(X1),r0));
        negative_grad = (gradsum < 0);
    }
    return bound;
}

//Brents Method for line-search using fixed number of iterations.
INLINE real_t BrentsMethod(coord3d* X, coord3d& r0, coord3d* X1, coord3d* X2){
    real_t a,b,s,d;
    a = (real_t)0.0; b = (real_t)1.0; 

    //To match python reference implementation by Buster.
    //b = FindLineSearchBound(X,r0,X1);

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
        /** Inverse quadratic interpolation **/
        if ( (f_a != f_c) && (f_b != f_c) )
        {
            s = a*f_a*f_c / ((f_a - f_b)*(f_a - f_c)) + b * f_a * f_c / ((f_b-f_a)*(f_b-f_c)) + c*f_a*f_b/((f_c-f_a)*(f_c-f_b));
        }else /** Secant Method **/
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
            s = (a+b) / (real_t)2.0; /** Bisection Method **/
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
INLINE real_t GSS(coord3d* X, coord3d& r0, coord3d* X1, coord3d* X2) const{
    constexpr real_t tau = (real_t)0.6180339887;
    //Line search x - values;
    real_t a = 0.0; real_t b = (real_t)1.0;
    
    //Use this bound to match reference python implementation by Buster.
    //b =FindLineSearchBound(X, r0, X1);
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

INLINE  void CG(coord3d* X, coord3d* X1, coord3d* X2, const size_t MaxIter)
{
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

    for (size_t i = 0; i < MaxIter; i++)
    {   
        alpha = GSS(X,s,X1,X2);
        if (alpha > (real_t)0.0){X1[node_id] = X[node_id] + alpha * s;}
        g1 = gradient(X1);
        //Polak Ribiere method
        g0_norm2 = reduction(sdata, dot(g0, g0));
        beta = max(reduction(sdata, dot(g1, (g1 - g0))) / g0_norm2,(real_t)0.0);

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

__global__ void kernel_optimize_batch(IsomerspaceForcefield::DeviceGraph G, const size_t max_iter){
    typedef device_coord3d coord3d;
    extern __shared__ device_real_t smem[];
    clear_cache(smem,Block_Size_Pow_2);
    
    device_real_t* base_pointer = smem + Block_Size_Pow_2;
    size_t offset               = blockIdx.x * blockDim.x;
    size_t node_id              = threadIdx.x;
    size_t N                    = blockDim.x;
    
    //Set VRAM pointer to start of each fullerene, as opposed to at the start of the isomerbatch.
    coord3d* X = reinterpret_cast<coord3d*>(G.X+3*offset);

    //Assign a section of L1 cache to each set of cartesian coordinates X, X1 and X2.
    coord3d* sX =reinterpret_cast<coord3d*>(base_pointer);
    coord3d* X1 =reinterpret_cast<coord3d*>(base_pointer+3*N);
    coord3d* X2 =reinterpret_cast<coord3d*>(base_pointer+6*N);  

                                
    sX[node_id] = X[node_id];   //Copy cartesian coordinates from DRAM to L1 Cache.
    X           = &sX[0];       //Switch coordinate pointer from DRAM to L1 Cache.

    //Pre-compute force constants and store in registers.
    Constants constants = Constants(G);
    NodeGraph nodeG     = NodeGraph(G);
    //NodeGraph bookkeeping = NodeGraph(&neighbours[0],&face_right[0],&next_on_face[0],&prev_on_face[0]);   

    //Create forcefield struct and use optimization algorithm to optimize the fullerene 
    ForceField FF = ForceField(nodeG, constants, smem);
    FF.CG(X,X1,X2,max_iter);
    BLOCK_SYNC;
    //Copy data back from L1 cache to VRAM 
    reinterpret_cast<coord3d*>(G.X)[offset + threadIdx.x]= X[threadIdx.x];   
}



__global__ void kernel_check_batch(IsomerspaceForcefield::DeviceGraph G, device_real_t* global_reduction_array){
    extern __shared__ device_real_t smem[];

    device_real_t* base_pointer = smem + Block_Size_Pow_2;
    clear_cache(smem,Block_Size_Pow_2);

    device_node_t node_id = threadIdx.x;
    size_t offset = blockIdx.x * blockDim.x;
    size_t N = blockDim.x;

    device_coord3d* X   = reinterpret_cast<device_coord3d*>(G.X+3*offset);
    device_coord3d* sX  = reinterpret_cast<device_coord3d*>(base_pointer);

    sX[node_id] = X[node_id];
    X           = &sX[0];

    Constants constants     = Constants(G);
    NodeGraph node_graph    = NodeGraph(G);
    ForceField FF           = ForceField(node_graph, constants, smem);

    device_real_t NodeEnergy = 0.0,         NodeBond_Error = 0.0,     NodeAngle_Error = 0.0,    NodeDihedral_Error = 0.0;
    device_real_t ArcEnergy  = 0.0,         ArcBond_Error = 0.0,      ArcAngle_Error = 0.0,     ArcDihedral_Error = 0.0;
    device_real_t NodeTotBondError = 0.0,   NodeTotAngleError = 0.0,  NodeTotDihedralError = 0.0;
    GRID_SYNC;

    for (uint8_t j = 0; j < 3; j++){
        auto arc            = ForceField::ArcData(j, X, node_graph);
        ArcEnergy           = arc.dihedral_energy(constants);
        ArcBond_Error       = abs(abs(arc.bond()      - d_get(constants.r0,j))          /d_get(constants.r0,j));
        ArcAngle_Error      = abs(abs(arc.angle()     - d_get(constants.angle0,j))      /d_get(constants.angle0,j));
        ArcDihedral_Error   = abs(abs(arc.dihedral()  - d_get(constants.inner_dih0,j))  /d_get(constants.inner_dih0,j));
        NodeEnergy += ArcEnergy; NodeBond_Error = max(ArcBond_Error,NodeBond_Error); NodeAngle_Error = max(ArcAngle_Error,NodeAngle_Error); NodeDihedral_Error = max(ArcDihedral_Error,NodeDihedral_Error);
        NodeTotAngleError += ArcBond_Error; NodeTotDihedralError += ArcDihedral_Error; NodeTotBondError += ArcBond_Error;
    }
    device_real_t Energy                = FF.energy(X);
    device_real_t GradNorm              = reduction(smem,dot(FF.gradient(X),FF.gradient(X)))/(device_real_t)blockDim.x;
    device_real_t MaxBond_Error         = reduction_max(reinterpret_cast<device_real_t*>(smem), NodeBond_Error);
    device_real_t MaxAngle_Error        = reduction_max(reinterpret_cast<device_real_t*>(smem), NodeAngle_Error);
    device_real_t MaxDihedral_Error     = reduction_max(reinterpret_cast<device_real_t*>(smem), NodeDihedral_Error);
    device_real_t RMS_Bond_Error        = sqrt(reduction(reinterpret_cast<device_real_t*>(smem),NodeTotBondError*NodeTotBondError)/blockDim.x);
    device_real_t RMS_Angle_Error       = sqrt(reduction(reinterpret_cast<device_real_t*>(smem),NodeTotAngleError*NodeTotAngleError)/blockDim.x);
    device_real_t RMS_Dihedral_Error    = sqrt(reduction(reinterpret_cast<device_real_t*>(smem),NodeTotDihedralError*NodeTotDihedralError)/blockDim.x);
    
    bool converged                      = (MaxBond_Error < 1e-1) && (MaxAngle_Error < 1.9e-1) && (MaxDihedral_Error < 1e-1) && (MaxBond_Error > 0.0) && (MaxAngle_Error > 0.0) && (MaxDihedral_Error > 0.0);
    global_reduction_array[blockIdx.x]  = (device_real_t)converged;
    
    for (size_t i = 0; i < gridDim.x; i++){
        GRID_SYNC;
        if(blockIdx.x == 0 && threadIdx.x == 0 && global_reduction_array[0]<0.5){printf("                    Error Summary                     \n====================================================\n Dihedral Max/RMS: \t%e | %e\n AngleMaxErr Max/RMS: \t%e | %e \n BondMaxErr Max/RMS: \t%e | %e \n \t\t   \t\t\t\t \n Energy/ Grad: \t\t%e | %e  \t\t \n====================================================\n\n", MaxDihedral_Error, RMS_Dihedral_Error, MaxAngle_Error, RMS_Angle_Error, MaxBond_Error, RMS_Bond_Error, Energy/blockDim.x, GradNorm);}
    }
    GRID_SYNC;
    device_real_t num_of_success = 0;
    if ((threadIdx.x + blockIdx.x * blockDim.x) == 0)
    {
        for (size_t i = 0; i < gridDim.x; i++)
        {
            num_of_success += global_reduction_array[i];
        }
    }
    GRID_SYNC;
    if(threadIdx.x + blockIdx.x == 0){printf("%d", (int)num_of_success); printf("/ %d Fullerenes Converged in Batch \n", (int)gridDim.x);}
}

__global__ void kernel_internal_coordinates(IsomerspaceForcefield::DeviceGraph G, IsomerspaceForcefield::InternalCoordinates c){
    size_t offset = blockIdx.x * blockDim.x;
    device_coord3d* X                       = &reinterpret_cast<device_coord3d*>(G.X)[offset];
    Constants constants     = Constants(G);
    NodeGraph node_graph    = NodeGraph(G);

    size_t tid = threadIdx.x + blockDim.x*blockIdx.x;
    for (uint8_t j = 0; j < 3; j++)
    {   
        auto arc                = ForceField::ArcData(j, X, node_graph);
        c.bonds[tid*3 + j]      = arc.bond();
        c.angles[tid*3 + j]     = arc.angle();
        c.dihedrals[tid*3 + j]  = arc.dihedral();
    }
}




void IsomerspaceForcefield::check_batch(){
    void* kernelArgs[] = {(void*)&d_graph, (void*)&global_reduction_array};
    cudaLaunchCooperativeKernel((void*)kernel_check_batch, dim3(batch_size,1,1), dim3(N,1,1), kernelArgs, sizeof(device_coord3d)*3*N + sizeof(device_real_t)*Block_Size_Pow_2);
    cudaDeviceSynchronize();
}

size_t IsomerspaceForcefield::get_batch_capacity(size_t N){
    cudaDeviceProp properties;
    cudaGetDeviceProperties(&properties,0);

    /** Compiling with --maxrregcount=64   is necessary to easily (singular blocks / fullerene) parallelize fullerenes of size 20-1024 !**/
    int fullerenes_per_SM;
    
    /** Needs 3 storage arrays for coordinates and 1 for reductions **/
    int sharedMemoryPerBlock = sizeof(device_coord3d)* 3 * N + sizeof(device_real_t)*Block_Size_Pow_2;

    /** Calculates maximum number of resident fullerenes on a single Streaming Multiprocessor, multiply with multi processor count to d_get total batch size**/
    cudaOccupancyMaxActiveBlocksPerMultiprocessor(&fullerenes_per_SM, kernel_optimize_batch, N, sharedMemoryPerBlock);

    return (size_t)(properties.multiProcessorCount*fullerenes_per_SM);
}

void IsomerspaceForcefield::get_cartesian_coordinates(device_real_t* X){
    cudaMemcpy(X, d_graph.X, sizeof(device_coord3d)*N*batch_size, cudaMemcpyDeviceToHost);
}

void IsomerspaceForcefield::get_internal_coordinates(device_real_t* bonds, device_real_t* angles, device_real_t* dihedrals){
    void* kernelArgs[] = {(void*)&d_graph, (void*)&d_coords};
    cudaLaunchCooperativeKernel((void*)kernel_internal_coordinates, dim3(batch_size,1,1), dim3(N,1,1), kernelArgs, sizeof(device_coord3d)*3*N + sizeof(device_real_t)*Block_Size_Pow_2);
    cudaMemcpy(bonds,       d_coords.bonds,     sizeof(device_coord3d)*N*batch_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(angles,      d_coords.angles,    sizeof(device_coord3d)*N*batch_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(dihedrals,   d_coords.dihedrals, sizeof(device_coord3d)*N*batch_size, cudaMemcpyDeviceToHost);
}

void IsomerspaceForcefield::optimize_batch(const size_t MaxIter){
    d_graph.copy_to_gpu(h_graph);

    printLastCudaError("Memcpy Failed! \n");
    auto start = std::chrono::system_clock::now();
    void* kernelArgs[] = {(void*)&d_graph,(void*)&MaxIter};
    
    cudaLaunchCooperativeKernel((void*)kernel_optimize_batch, dim3(batch_size, 1, 1), dim3(N, 1, 1), kernelArgs, sizeof(device_coord3d)* 3 * N + sizeof(device_real_t)*Block_Size_Pow_2);
    cudaDeviceSynchronize();
    auto end = std::chrono::system_clock::now();
    printLastCudaError("Kernel launch failed: ");
    std::cout << "Elapsed time: " << (end-start)/ 1ms << "ms\n" ;
    std::cout << "Estimated Performance " << ((device_real_t)(batch_size)/(std::chrono::duration_cast<std::chrono::microseconds>(end-start)).count()) * 1.0e6 << "Fullerenes/s \n";
}

void IsomerspaceForcefield::insert_isomer_batch(const DeviceGraph& G){
    h_graph = G;
    batch_size += G.batch_size;
}

/*
void IsomerspaceForcefield::insert_isomer(const FullereneGraph& G, const vector<coord3d> &X0){
    size_t offset = batch_size*3*N;
    for (device_node_t u = 0; u < N; u++){
        for (int j = 0; j < 3; j++){
            device_node_t v = G.neighbours[u][j];
            size_t arc_index = u*3 + j + offset;
            h_graph.neighbours  [arc_index] = v;
            h_graph.next_on_face[arc_index] = G.next_on_face(u,v);
            h_graph.prev_on_face[arc_index] = G.prev_on_face(u,v);
            h_graph.face_right  [arc_index] = G.face_size(u,v);
            h_graph.X           [arc_index] = X0[u][j];
        }   
    }   
    batch_size++;
}*/

IsomerspaceForcefield::IsomerspaceForcefield(const size_t N)
{
    this->batch_capacity = get_batch_capacity((size_t)N);
    std::cout << "\nIsomerspace Capacity: " << this->batch_capacity << "\n";
    this->N = N;
    d_coords.allocate(N,batch_capacity);
    d_graph.allocate(N,batch_capacity);
    h_graph.allocate_host(N,batch_capacity);
    cudaMalloc(&global_reduction_array, sizeof(device_real_t)*N*batch_capacity);
    printLastCudaError("Kernel class instansiation failed!");
}

IsomerspaceForcefield::~IsomerspaceForcefield()
{   
    //Destructors of member objects ensures that cudaFree and free is called respectively on all allocated arrays.
    d_graph.free();
    d_coords.free();
    h_graph.free_host();
    cudaDeviceReset();
}



