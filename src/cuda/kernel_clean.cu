#include "fullerenes/gpu/isomerspace_forcefield.hh"
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

namespace IsomerspaceForcefield {

typedef device_real_t real_t;
typedef device_node_t node_t;

#include "coord3d.cu"
#include "helper_functions.cu"


using namespace std::literals;
namespace cg = cooperative_groups;

struct DevicePointers; struct HostPointers; struct CoordinatePointers; struct HarmonicConstantPointers;

/** This struct was made to reduce signature cluttering of device functions, it is simply a container for default arguments which are shared between functions**/
template <int Block_Size_Pow_2>
struct ForceField{
    coord3d r1,              //d_i+1 In conjugated gradient algorithm  Buster p. 36
            r0;                  //d_0 In conjugated gradient algorithm Buster p. 36

    const BookkeepingData bdat;         //Contains face-information and neighbour-information. Both of which are constant in the lifespan of this struct. 
    const Constants<coord3d> constants; //Contains force-constants and equillibrium-parameters. Constant in the lifespan of this struct.

    cg::thread_group group_handle = cg::this_thread_block(); //Synchronization handle may be changed to grid-wide synchronization 
    node_t node_id; //In case the code needs to be turned into grid-sized fullerenes. Parallel lockstep implementation of this is nontrivial to optimize in terms of blocksizes.
    real_t* sdata;  //Pointer to start of L1 cache array, used exclusively for reduction.

    __device__ ForceField(  const BookkeepingData &b,
                            const Constants<coord3d> &c, 
                            real_t* sdata,
                            node_t node_id): bdat(b), constants(c), sdata(sdata), node_id(node_id) {}


//Container for all energy and gradient evaluations with respect to an arc, eg. AB, AC or AD.
struct ArcData{
    //124 FLOPs;

    __device__ ArcData(const node_t a, const uint8_t j, const coord3d* __restrict__ X, const BookkeepingData& bdat){   
        this->j = j;   

        real_t r_rmp;
        coord3d ap, am, ab, ac, ad, mp;
        coord3d X_a = X[a]; coord3d X_b = X[bdat.neighbours[j]];
        //printf("Index: %d \n", a*3 + j);

        //Compute the arcs ab, ac, ad, bp, bm, ap, am, mp, bc and cd
        ab = (X_b - X_a);  r_rab = bond_length(ab); ab_hat = r_rab * ab;
        ac = (X[bdat.neighbours[(j+1)%3]] - X_a); r_rac = bond_length(ac); ac_hat = r_rac * ac; rab = non_resciprocal_bond_length(ab);
        ad = (X[bdat.neighbours[(j+2)%3]] - X_a); r_rad = bond_length(ad); ad_hat = r_rad * ad;
        
        coord3d bp = (X[bdat.next_on_face[j]] - X_b); bp_hat = unit_vector(bp);
        coord3d bm = (X[bdat.prev_on_face[j]] - X_b); bm_hat = unit_vector(bm);

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
    INLINE coord3d inner_angle_gradient(const Constants<coord3d>& c) const
    {   
        real_t cos_angle = angle(); //Inner angle of arcs ab,ac.
        coord3d grad = cos_angle * (ab_hat * r_rab + ac_hat * r_rac) - ab_hat * r_rac - ac_hat* r_rab; //Derivative of inner angle: Eq. 21. 
        return get(c.f_inner_angle,j) * harmonic_energy_gradient(get(c.angle0,j), cos_angle, grad); //Harmonic Energy Gradient: Eq. 21. multiplied by harmonic term.
    }
    //Computes gradient related to bending of outer angles. ~20 FLOPs
    INLINE coord3d outer_angle_gradient_m(const Constants<coord3d>& c) const
    {
        real_t cos_angle = -dot(ab_hat, bm_hat); //Compute outer angle. ab,bm
        coord3d grad = (bm_hat + ab_hat * cos_angle) * r_rab; //Derivative of outer angles Eq. 30. Buster Thesis
        return get(c.f_outer_angle_m,j) * harmonic_energy_gradient(get(c.outer_angle_m0,j),cos_angle,grad); //Harmonic Energy Gradient: Eq. 30 multiplied by harmonic term.
    }
    INLINE coord3d outer_angle_gradient_p(const Constants<coord3d>& c) const
    {
        real_t cos_angle = -dot(ab_hat, bp_hat); //Compute outer angle. ab,bp
        coord3d grad = (bp_hat + ab_hat * cos_angle) * r_rab; //Derivative of outer angles Eq. 28. Buster Thesis
        return get(c.f_outer_angle_p,j) * harmonic_energy_gradient(get(c.outer_angle_p0,j),cos_angle,grad); //Harmonic Energy Gradient: Eq. 28 multiplied by harmonic term.
    }
    // Chain rule terms for dihedral calculation
    //Computes gradient related to dihedral/out-of-plane term. ~75 FLOPs
    INLINE coord3d inner_dihedral_gradient(const Constants<coord3d>& c) const
    {
        coord3d nabc, nbcd; real_t cos_b, cos_c, r_sin_b, r_sin_c;
        cos_b = dot(ba_hat,bc_hat); r_sin_b = (real_t)1.0/sqrt((real_t)1.0 - cos_b*cos_b); nabc = cross(ba_hat, bc_hat) * r_sin_b;
        cos_c = dot(-bc_hat,cd_hat); r_sin_c = (real_t)1.0/sqrt((real_t)1.0 - cos_c*cos_c); nbcd = cross(-bc_hat,cd_hat) * r_sin_c;

        real_t cos_beta = dot(nabc, nbcd); //Inner dihedral angle from planes abc,bcd.
        real_t cot_b = cos_b * r_sin_b * r_sin_b; //cos(b)/sin(b)^2

        //Derivative w.r.t. inner dihedral angle F and G in Eq. 26
        coord3d grad = cross(bc_hat, nbcd) * r_sin_b * r_rab - ba_hat * cos_beta * r_rab + (cot_b * cos_beta * r_rab) * (bc_hat - ba_hat * cos_b);

        return get(c.f_inner_dihedral,j) * harmonic_energy_gradient(get(c.inner_dih0,j), cos_beta, grad); //Eq. 26.
    }

    

    //Computes gradient from dihedral angles constituted by the planes bam, amp ~162 FLOPs
    INLINE coord3d outer_dihedral_gradient_a(const Constants<coord3d>& c) const
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

        return get(c.f_outer_dihedral,j) * harmonic_energy_gradient(get(c.outer_dih0_a,j), cos_beta, grad);
    }

    //Computes gradient from dihedral angles constituted by the planes nbmp, nmpa ~92 FLOPs
    INLINE coord3d outer_dihedral_gradient_m(const Constants<coord3d>& c) const
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
        return get(c.f_outer_dihedral,j) * harmonic_energy_gradient(get(c.outer_dih0_m,j), cos_beta, grad);
    }

    //Computes gradient from dihedral angles constituted by the planes bpa, pam ~162 FLOPs
    INLINE coord3d outer_dihedral_gradient_p(const Constants<coord3d>& c) const
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
        return get(c.f_outer_dihedral,j) * harmonic_energy_gradient(get(c.outer_dih0_p,j), cos_beta, grad);
    }
    // Internal coordinate gradients
    INLINE coord3d bond_length_gradient(const Constants<coord3d>& c) const { return get(c.f_bond,j) * harmonic_energy_gradient(bond(),get(c.r0,j),ab_hat);}
    //Sum of angular gradient components.
    INLINE coord3d angle_gradient(const Constants<coord3d>& c) const { return inner_angle_gradient(c) + outer_angle_gradient_p(c) + outer_angle_gradient_m(c);}
    //Sum of inner and outer dihedral gradient components.
    INLINE coord3d dihedral_gradient(const Constants<coord3d>& c) const { return inner_dihedral_gradient(c) + outer_dihedral_gradient_a(c) + outer_dihedral_gradient_m(c) + outer_dihedral_gradient_p(c);}
    //coord3d flatness()             const { return ;  }   
    
    INLINE real_t bond_energy(const Constants<coord3d>& c) const {return (real_t)0.5 *get(c.f_bond,j) *harmonic_energy(bond(),get(c.r0,j));}
    INLINE real_t bend_energy(const Constants<coord3d>& c) const {return get(c.f_inner_angle,j)* harmonic_energy(angle(),get(c.angle0,j));}
    INLINE real_t dihedral_energy(const Constants<coord3d>& c) const {return get(c.f_inner_dihedral,j)* harmonic_energy(dihedral(),get(c.inner_dih0,j));}
    //Harmonic energy contribution from bond stretching, angular bending and dihedral angle bending.
    //71 FLOPs
    INLINE real_t energy(const Constants<coord3d>& c) const {return bond_energy(c) + bend_energy(c) + dihedral_energy(c); }
    //Sum of bond, angular and dihedral gradient components.
    INLINE coord3d gradient(const Constants<coord3d>& c) const{return bond_length_gradient(c) + angle_gradient(c) + dihedral_gradient(c);}

    
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
    cg::sync(group_handle);
    coord3d grad = {0.0, 0.0, 0.0};

    for (uint8_t j = 0; j < 3; j++ ){
        ArcData arc = ArcData(node_id, j, X, bdat);
        grad += arc.gradient(constants);
        //set(grad,j,get(constants.f_outer_dihedral,j));
    }

    return grad;
}

INLINE real_t energy(coord3d* X) const {
    cg::sync(group_handle);
    real_t arc_energy = (real_t)0.0;

    //(71 + 124) * 3 * N  = 585*N FLOPs
    for (uint8_t j = 0; j < 3; j++ ){
        ArcData arc = ArcData(node_id, j, X, bdat);
        arc_energy += arc.energy(constants);
    }
    
    return reduction<Block_Size_Pow_2>(sdata, arc_energy);;

}
INLINE real_t FindLineSearchBound(coord3d* X, coord3d& r0, coord3d* X1){
    real_t bound = 1e-5;
    bool negative_grad = true;
    while (negative_grad)
    {   
        bound *= (real_t)1.5;
        X1[node_id] = X[node_id] + bound * r0;
        real_t gradsum = reduction<Block_Size_Pow_2>(sdata, dot(gradient(X1),r0));
        negative_grad = (gradsum < 0);
    }
    return bound;
}

INLINE real_t BrentsMethod(coord3d* X, coord3d& r0, coord3d* X1, coord3d* X2){
    real_t a,b,s,d;
    a = (real_t)0.0; b = (real_t)1.0; 

    //To match python reference implementation by Buster.
    //b = FindLineSearchBound(X,r0,X1);

    X1[node_id] = X[node_id] + a * r0;
    X2[node_id] = X[node_id] + b * r0;

    real_t f_a = energy(X1);
    real_t f_b = energy(X2);

    if (f_a < f_b){swap(a,b);}

    real_t c = a; real_t f_c = f_a;
    bool flag = true;
    int Iterations = 0;

    while (abs(b-a) > (real_t)5e-8)
    {   
        Iterations++;

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
        real_t f_s = energy(X1);
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
            swap(a,b); swap(f_a,f_b);
        }
    }
    real_t alpha = s;

    return alpha;
}

INLINE real_t GSS(coord3d* X, coord3d& r0, coord3d* X1, coord3d* X2){
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

    for (uint8_t i = 0; i < 30; i++){
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
    //if (f1 > energy(X)) {return (real_t)0.0;}
    //Line search coefficient
    real_t alpha = (a+b)/(real_t)2.0;
    return alpha;
}

INLINE  void CG(coord3d* X, coord3d* X1, coord3d* X2, const size_t MaxIter)
{
    real_t alpha, beta, g0_norm2, s_norm;
    coord3d g0,g1,s;
    g0 = -gradient(X); 
    s = g0;

    //Normalize To match reference python implementation by Buster.
    s_norm = sqrt(reduction<Block_Size_Pow_2>(sdata, dot(s,s)));
    s /= s_norm;

    for (size_t i = 0; i < MaxIter; i++)
    {   
        //print_single(s_norm);
        alpha = GSS(X,s,X1,X2);
        //print_single(s_norm);
        //if (alpha > (real_t)0.0){X1[node_id] = X[node_id] + alpha * s;}
        X1[node_id] = X[node_id] + alpha * s;
        g1 = gradient(X1);
        //Polak Ribiere method
        g0_norm2 = reduction<Block_Size_Pow_2>(sdata, dot(g0, g0));
        beta = max(reduction<Block_Size_Pow_2>(sdata, dot(g1, (g1 - g0))) / g0_norm2,(real_t)0.0);
        if (energy(X1) > energy(X)){X1[node_id] = X[node_id]; g1 = g0; beta = (real_t) 0.0;}
        //if (alpha > (real_t)0.0){X[node_id] = X1[node_id]; g1 = g0; beta = (real_t) 0.0;}

        s = -g1 + beta*s;
        g0 = g1;
        X[node_id] = X1[node_id];

        //Normalize Search Direction using MaxNorm or 2Norm
        //s_norm = reduction_max<Block_Size_Pow_2>(sdata, max(max(s.x,s.y),s.z));
        s_norm = sqrt(reduction<Block_Size_Pow_2>(sdata, dot(s,s)));
        s /= s_norm;
    }   
} 


};

template <int Block_Size_Pow_2>
__global__ void kernel_OptimizeBatch(DevicePointers p, const size_t N, const size_t MaxIter){
    extern __shared__ real_t smem[];
    real_t* base_pointer = smem + Block_Size_Pow_2;
    cg::grid_group grid = cg::this_grid();
    cg::thread_block block = cg::this_thread_block();

    node_t node_id = threadIdx.x;
    size_t offset = blockIdx.x * blockDim.x;

    coord3d* X = reinterpret_cast<coord3d*>(p.X+3*offset);
    coord3d* X1 = reinterpret_cast<coord3d*>(p.X1+3*offset);
    coord3d* X2 = reinterpret_cast<coord3d*>(p.X2+3*offset);
    
    coord3d* sX = reinterpret_cast<coord3d*>(base_pointer);
    coord3d* sX1 =reinterpret_cast<coord3d*>(base_pointer+3*N);
    coord3d* sX2 =reinterpret_cast<coord3d*>(base_pointer+6*N);  

    cg::sync(grid);

    sX[node_id] = X[node_id];
    sX1[node_id] = sX[node_id];
    X = &sX[0]; X1 = &sX1[0]; X2 = &sX2[0];
    
    //Pre-compute force constants and store in registers.
    BookkeepingData bookit = BookkeepingData(&p.neighbours[3*offset],&p.face_right[3*offset],&p.next_on_face[3*offset],&p.prev_on_face[3*offset]);
    Constants<coord3d> constants = compute_constants<coord3d>(bookit, node_id);

    //Load constant bookkeeping data into registers.
    const node_t neighbours[3] = {p.neighbours[3*(offset+node_id)],p.neighbours[3*(offset+node_id) + 1],p.neighbours[3*(offset+node_id) + 2]};
    const uint8_t face_right[3] = {p.face_right[3*(offset+node_id)],p.face_right[3*(offset+node_id) + 1],p.face_right[3*(offset+node_id) + 2]};;
    const node_t next_on_face[3] = {p.next_on_face[3*(offset+node_id)],p.next_on_face[3*(offset+node_id) + 1],p.next_on_face[3*(offset+node_id) + 2]};
    const node_t prev_on_face[3] = {p.prev_on_face[3*(offset+node_id)],p.prev_on_face[3*(offset+node_id) + 1],p.prev_on_face[3*(offset+node_id) + 2]};
    BookkeepingData bookkeeping = BookkeepingData(&neighbours[0],&face_right[0],&next_on_face[0],&prev_on_face[0]);   



    cg::sync(grid);

    /** Create forcefield struct and use optimization algorithm to optimize the fullerene **/
    ForceField<Block_Size_Pow_2> FF = ForceField<Block_Size_Pow_2>(bookkeeping, constants, smem, node_id);
    FF.CG(X,X1,X2,MaxIter);
    cg::sync(grid);
    /** Copy data back from L1 cache to VRAM **/    
    reinterpret_cast<coord3d*>(p.X)[offset + threadIdx.x] = X[threadIdx.x];
}

template <int Block_Size_Pow_2>
__global__ void FullereneProperties(DevicePointers p){
    extern __shared__ coord3d sdata[];
    size_t offset = blockIdx.x * blockDim.x;
    coord3d* X = &reinterpret_cast<coord3d*>(p.X)[offset];


    BookkeepingData bookit = BookkeepingData(&p.neighbours[3*(offset)],&p.face_right[3*(offset)],&p.next_on_face[3*(offset)],&p.prev_on_face[3*(offset)]);
    Constants<coord3d> constants = compute_constants<coord3d>(bookit,threadIdx.x);
    BookkeepingData bdat = BookkeepingData(&p.neighbours[3*(offset + threadIdx.x)],&p.face_right[3*(offset + threadIdx.x)],&p.next_on_face[3*(offset + threadIdx.x)],&p.prev_on_face[3*(offset + threadIdx.x)]);

    //Use X1 buffer as storage array.
    coord3d* NodeEnergyCoord = &reinterpret_cast<coord3d*>(p.X1)[offset];
    real_t NodeEnergy = (real_t)0.0; real_t NodeBond_Error = (real_t)0.0; real_t NodeAngle_Error = (real_t)0.0; real_t NodeDihedral_Error = (real_t)0.0;
    real_t ArcEnergy = (real_t)0.0; real_t ArcBond_Error = (real_t)0.0; real_t ArcAngle_Error = (real_t)0.0; real_t ArcDihedral_Error = (real_t)0.0;
    real_t NodeTotBondError = 0.0, NodeTotAngleError= 0.0, NodeTotDihedralError = 0.0;
    ForceField<Block_Size_Pow_2> FF = ForceField<Block_Size_Pow_2>(bdat, constants, reinterpret_cast<real_t*>(sdata), threadIdx.x);
    
    for (uint8_t j = 0; j < 3; j++)
    {
        auto arc = ForceField<Block_Size_Pow_2>::ArcData(threadIdx.x, j, X, bdat);
        ArcEnergy = arc.dihedral_energy(constants);

        ArcBond_Error = abs(abs(arc.bond() - get(constants.r0,j))/get(constants.r0,j));
        ArcAngle_Error =  abs(abs(arc.angle() - get(constants.angle0,j))/get(constants.angle0,j));
        ArcDihedral_Error = abs(abs(arc.dihedral() - get(constants.inner_dih0,j))/get(constants.inner_dih0,j));

        NodeEnergy += ArcEnergy; NodeBond_Error = max(ArcBond_Error,NodeBond_Error); NodeAngle_Error = max(ArcAngle_Error,NodeAngle_Error); NodeDihedral_Error = max(ArcDihedral_Error,NodeDihedral_Error);
        NodeTotAngleError += ArcBond_Error; NodeTotDihedralError += ArcDihedral_Error; NodeTotBondError += ArcBond_Error;

    }
    real_t Energy = FF.energy(X);
    real_t MaxBond_Error =  reduction_max<Block_Size_Pow_2>(reinterpret_cast<real_t*>(sdata), NodeBond_Error);
    real_t MaxAngle_Error =  reduction_max<Block_Size_Pow_2>(reinterpret_cast<real_t*>(sdata), NodeAngle_Error);
    real_t MaxDihedral_Error =  reduction_max<Block_Size_Pow_2>(reinterpret_cast<real_t*>(sdata), NodeDihedral_Error);
    
    real_t RMS_Bond_Error = sqrt(reduction<Block_Size_Pow_2>(reinterpret_cast<real_t*>(sdata), NodeTotBondError*NodeTotBondError)/blockDim.x);
    real_t RMS_Angle_Error = sqrt(reduction<Block_Size_Pow_2>(reinterpret_cast<real_t*>(sdata), NodeTotAngleError*NodeTotAngleError)/blockDim.x);
    real_t RMS_Dihedral_Error = sqrt(reduction<Block_Size_Pow_2>(reinterpret_cast<real_t*>(sdata), NodeTotDihedralError*NodeTotDihedralError)/blockDim.x);
    
    
    if (threadIdx.x == 0){

    if ((MaxBond_Error < 1e-1) && (MaxAngle_Error < 1e-1) && (MaxDihedral_Error < 1e-1))
    {   
        p.gdata[blockIdx.x] = (real_t)1.0;

    } else
    {
        p.gdata[blockIdx.x] = (real_t)0.0;
    }}
    cg::sync(cg::this_grid());
    if(threadIdx.x + blockIdx.x == 0){printf("                    Error Summary                     \n====================================================\n Dihedral Max/RMS: \t%e | %e\n AngleMaxErr Max/RMS: \t%e | %e \n BondMaxErr Max/RMS: \t%e | %e \n \t\t   \t\t\t\t \n Energy/mol: \t\t%e   \t\t \n====================================================\n\n", MaxDihedral_Error, RMS_Dihedral_Error, MaxAngle_Error, RMS_Angle_Error, MaxBond_Error, RMS_Bond_Error, Energy/blockDim.x);}



    cg::sync(cg::this_grid());
    real_t Success = 0;
    if ((threadIdx.x + blockIdx.x * blockDim.x) == 0)
    {
        for (size_t i = 0; i < gridDim.x; i++)
        {
            Success += p.gdata[i];
        }
    }
    cg::sync(cg::this_grid());
    if(threadIdx.x + blockIdx.x == 0){printf("%d", (int)Success); printf("/ %d Fullerenes Converged in Batch \n", (int)gridDim.x);}

}

template <int Block_Size_Pow_2>
__global__ void kernel_InternalCoordinates(DevicePointers p, CoordinatePointers c){
    size_t offset = blockIdx.x * blockDim.x;
    coord3d* X = &reinterpret_cast<coord3d*>(p.X)[offset];
    BookkeepingData bookit = BookkeepingData(&p.neighbours[3*(offset)],&p.face_right[3*(offset)],&p.next_on_face[3*(offset)],&p.prev_on_face[3*(offset)]);
    Constants<coord3d> constants = compute_constants<coord3d>(bookit,threadIdx.x);
    BookkeepingData bdat = BookkeepingData(&p.neighbours[3*(offset + threadIdx.x)],&p.face_right[3*(offset + threadIdx.x)],&p.next_on_face[3*(offset + threadIdx.x)],&p.prev_on_face[3*(offset + threadIdx.x)]);

    size_t tid = threadIdx.x + blockDim.x*blockIdx.x;
    for (uint8_t j = 0; j < 3; j++)
    {   
        auto arc = ForceField<Block_Size_Pow_2>::ArcData(threadIdx.x, j, X, bdat);
        c.bonds[tid*3 + j] = arc.bond();
        c.angles[tid*3 + j] = arc.angle();
        c.outer_angles_m[tid*3 + j] = arc.outer_angle_m();
        c.outer_angles_p[tid*3 + j] = arc.outer_angle_p();
        c.dihedrals[tid*3 + j] = arc.dihedral();
        c.outer_dihedrals_a[tid*3 + j] = arc.outer_dihedral_a();
        c.outer_dihedrals_m[tid*3 + j] = arc.outer_dihedral_m();
        c.outer_dihedrals_p[tid*3 + j] = arc.outer_dihedral_p();
    }
}
template <int Block_Size_Pow_2>
__global__ void kernel_HarmonicConstants(DevicePointers p, HarmonicConstantPointers c0){
    size_t offset = blockIdx.x * blockDim.x;
    BookkeepingData bookit = BookkeepingData(&p.neighbours[3*(offset)],&p.face_right[3*(offset)],&p.next_on_face[3*(offset)],&p.prev_on_face[3*(offset)]);
    Constants<coord3d> constants = compute_constants<coord3d>(bookit,threadIdx.x);
    BookkeepingData bdat = BookkeepingData(&p.neighbours[3*(offset + threadIdx.x)],&p.face_right[3*(offset + threadIdx.x)],&p.next_on_face[3*(offset + threadIdx.x)],&p.prev_on_face[3*(offset + threadIdx.x)]);

    size_t tid = threadIdx.x + blockDim.x*blockIdx.x;
    for (uint8_t j = 0; j < 3; j++)
    {   
        c0.bonds_0[tid*3 + j] = get(constants.r0,j);
        c0.angles_0[tid*3 + j] = get(constants.angle0,j);
        c0.outer_angles_m0[tid*3 + j] = get(constants.outer_angle_m0,j);
        c0.outer_angles_p0[tid*3 + j] = get(constants.outer_angle_p0,j);
        c0.dihedrals_0[tid*3 + j] = get(constants.inner_dih0,j);
        c0.outer_dihedrals_a0[tid*3 + j] = get(constants.outer_dih0_a,j);
        c0.outer_dihedrals_m0[tid*3 + j] = get(constants.outer_dih0_m,j);
        c0.outer_dihedrals_p0[tid*3 + j] = get(constants.outer_dih0_p,j);
    }
}

template <int Block_Size_Pow_2>
__global__ void kernel_Gradients(DevicePointers p){
    extern __shared__ real_t sdat[];
    size_t offset = blockIdx.x * blockDim.x;
    coord3d* X = &reinterpret_cast<coord3d*>(p.X)[offset];
    BookkeepingData bookit = BookkeepingData(&p.neighbours[3*(offset)],&p.face_right[3*(offset)],&p.next_on_face[3*(offset)],&p.prev_on_face[3*(offset)]);
    Constants<coord3d> constants = compute_constants<coord3d>(bookit,threadIdx.x);
    BookkeepingData bdat = BookkeepingData(&p.neighbours[3*(offset + threadIdx.x)],&p.face_right[3*(offset + threadIdx.x)],&p.next_on_face[3*(offset + threadIdx.x)],&p.prev_on_face[3*(offset + threadIdx.x)]);

    ForceField<Block_Size_Pow_2> FF = ForceField<Block_Size_Pow_2>(bdat, constants, sdat, threadIdx.x);
    size_t tid = threadIdx.x + blockDim.x*blockIdx.x;
    coord3d grad = FF.gradient(X);
    for (uint8_t j = 0; j < 3; j++)
    {   
        //auto arc = ForceField<Block_Size_Pow_2>::ArcData(threadIdx.x, j, X, bdat);
        //grad += arc.gradient(constants);
        //set(grad,j,get(constants.f_inner_dihedral,j));
    }
    reinterpret_cast<coord3d*>(p.gradients)[tid] = grad;
    real_t norm = reduction<Block_Size_Pow_2>(sdat,dot(grad,grad));
    print_single(sqrt(norm));


}

template <int Block_Size_Pow_2>    
size_t computeBatchSize(size_t N){
    cudaDeviceProp properties;
    cudaGetDeviceProperties(&properties,0);

    /** Compiling with --maxrregcount=64   is necessary to easily (singular blocks / fullerene) parallelize fullerenes of size 20-1024 !**/
    int fullerenes_per_SM;
    
    /** Needs 3 storage arrays for coordinates and 1 for reductions **/
    int sharedMemoryPerBlock = sizeof(coord3d)* 3 * (Block_Size_Pow_2 + 1) + sizeof(real_t)*N;

    /** Calculates maximum number of resident fullerenes on a single Streaming Multiprocessor, multiply with multi processor count to get total batch size**/
    cudaOccupancyMaxActiveBlocksPerMultiprocessor(&fullerenes_per_SM, kernel_OptimizeBatch<Block_Size_Pow_2>, N, sharedMemoryPerBlock);

    return (size_t)(properties.multiProcessorCount*fullerenes_per_SM);
}

void AllocateDevicePointers(DevicePointers& p, CoordinatePointers& c, HarmonicConstantPointers& c0, size_t N, size_t batch_size){
    cudaMalloc(&p.X, sizeof(coord3d)*N*batch_size);
    cudaMalloc(&p.X1, sizeof(coord3d)*N*batch_size);
    cudaMalloc(&p.X2, sizeof(coord3d)*N*batch_size);
    cudaMalloc(&p.neighbours, sizeof(node_t)*3*N*batch_size);
    cudaMalloc(&p.next_on_face, sizeof(node_t)*3*N*batch_size);
    cudaMalloc(&p.prev_on_face, sizeof(node_t)*3*N*batch_size);
    cudaMalloc(&p.face_right, sizeof(uint8_t)*3*N*batch_size);
    cudaMalloc(&p.gdata, sizeof(real_t)*batch_size);
    cudaMalloc(&p.gradients, sizeof(real_t)*3*N*batch_size);

    cudaMalloc(&c.bonds, sizeof(real_t)*3*N*batch_size);
    cudaMalloc(&c.angles, sizeof(real_t)*3*N*batch_size);
    cudaMalloc(&c.outer_angles_m, sizeof(real_t)*3*N*batch_size);
    cudaMalloc(&c.outer_angles_p, sizeof(real_t)*3*N*batch_size);
    cudaMalloc(&c.dihedrals, sizeof(real_t)*3*N*batch_size);
    cudaMalloc(&c.outer_dihedrals_a, sizeof(real_t)*3*N*batch_size);
    cudaMalloc(&c.outer_dihedrals_m, sizeof(real_t)*3*N*batch_size);
    cudaMalloc(&c.outer_dihedrals_p, sizeof(real_t)*3*N*batch_size);

    cudaMalloc(&c0.bonds_0, sizeof(real_t)*3*N*batch_size);
    cudaMalloc(&c0.angles_0, sizeof(real_t)*3*N*batch_size);
    cudaMalloc(&c0.outer_angles_m0, sizeof(real_t)*3*N*batch_size);
    cudaMalloc(&c0.outer_angles_p0, sizeof(real_t)*3*N*batch_size);
    cudaMalloc(&c0.dihedrals_0, sizeof(real_t)*3*N*batch_size);
    cudaMalloc(&c0.outer_dihedrals_a0, sizeof(real_t)*3*N*batch_size);
    cudaMalloc(&c0.outer_dihedrals_m0, sizeof(real_t)*3*N*batch_size);
    cudaMalloc(&c0.outer_dihedrals_p0, sizeof(real_t)*3*N*batch_size);
}

void FreePointers(DevicePointers& p, CoordinatePointers& c, HarmonicConstantPointers& c0){
    cudaFree(p.X);
    cudaFree(p.X1);
    cudaFree(p.X2);
    cudaFree(p.neighbours);
    cudaFree(p.next_on_face);
    cudaFree(p.prev_on_face);
    cudaFree(p.face_right);
    cudaFree(p.gdata);
    cudaFree(p.gradients);

    cudaFree(c.bonds);
    cudaFree(c.angles);
    cudaFree(c.outer_angles_m);
    cudaFree(c.outer_angles_p);
    cudaFree(c.dihedrals);
    cudaFree(c.outer_dihedrals_a);
    cudaFree(c.outer_dihedrals_m);
    cudaFree(c.outer_dihedrals_p);
    
    cudaFree(c0.bonds_0);
    cudaFree(c0.angles_0);
    cudaFree(c0.outer_angles_m0);
    cudaFree(c0.outer_angles_p0);
    cudaFree(c0.dihedrals_0);
    cudaFree(c0.outer_dihedrals_a0);
    cudaFree(c0.outer_dihedrals_m0);
    cudaFree(c0.outer_dihedrals_p0);
    
    cudaDeviceReset();
}

void CopyToDevice(DevicePointers& p, const HostPointers& h, const size_t N, const size_t batch_size){
    cudaMemcpy(p.X, h.h_X, sizeof(coord3d)*N*batch_size , cudaMemcpyHostToDevice);
    cudaMemcpy(p.neighbours, h.h_cubic_neighbours, sizeof(node_t)*3*N*batch_size, cudaMemcpyHostToDevice);
    cudaMemcpy(p.next_on_face, h.h_next_on_face, sizeof(node_t)*3*N*batch_size, cudaMemcpyHostToDevice);
    cudaMemcpy(p.prev_on_face, h.h_prev_on_face, sizeof(node_t)*3*N*batch_size, cudaMemcpyHostToDevice);
    cudaMemcpy(p.face_right, h.h_face_right, sizeof(uint8_t)*3*N*batch_size, cudaMemcpyHostToDevice);
}

template <int Block_Size_Pow_2>
void CheckBatch(DevicePointers& p, const HostPointers& h, const size_t N, const size_t batch_size){
    CopyToDevice(p,h,N,batch_size);
    void* kernelArgs[] = {(void*)&p};
    cudaLaunchCooperativeKernel((void*)FullereneProperties<Block_Size_Pow_2>, dim3(batch_size,1,1), dim3(N,1,1), kernelArgs, sizeof(coord3d)*3*N, NULL);
}

template <int Block_Size_Pow_2>
void InternalCoordinates(DevicePointers& p,const HostPointers& h, CoordinatePointers& h_c, CoordinatePointers& d_c, const size_t N, const size_t batch_size, real_t* bonds, real_t* angles, real_t* dihedrals){
    CopyToDevice(p,h,N,batch_size);
    void* kernelArgs[] = {(void*)&p, (void*)&d_c};
    cudaLaunchCooperativeKernel((void*)kernel_InternalCoordinates<Block_Size_Pow_2>, dim3(batch_size,1,1), dim3(N,1,1), kernelArgs, sizeof(coord3d)*3*N, NULL);
    cudaDeviceSynchronize();
    cudaMemcpy(h_c.bonds, d_c.bonds, sizeof(coord3d)*N*batch_size , cudaMemcpyDeviceToHost);
    cudaMemcpy(h_c.angles, d_c.angles, sizeof(coord3d)*N*batch_size , cudaMemcpyDeviceToHost);
    cudaMemcpy(h_c.outer_angles_m, d_c.outer_angles_m, sizeof(coord3d)*N*batch_size , cudaMemcpyDeviceToHost);
    cudaMemcpy(h_c.outer_angles_p, d_c.outer_angles_p, sizeof(coord3d)*N*batch_size , cudaMemcpyDeviceToHost);
    cudaMemcpy(h_c.dihedrals, d_c.dihedrals, sizeof(coord3d)*N*batch_size , cudaMemcpyDeviceToHost);
    cudaMemcpy(h_c.outer_dihedrals_a, d_c.outer_dihedrals_a, sizeof(coord3d)*N*batch_size , cudaMemcpyDeviceToHost);
    cudaMemcpy(h_c.outer_dihedrals_m, d_c.outer_dihedrals_m, sizeof(coord3d)*N*batch_size , cudaMemcpyDeviceToHost);
    cudaMemcpy(h_c.outer_dihedrals_p, d_c.outer_dihedrals_p, sizeof(coord3d)*N*batch_size , cudaMemcpyDeviceToHost);
}

template <int Block_Size_Pow_2>
void HarmonicConstants(DevicePointers& p, const HostPointers& h, HarmonicConstantPointers& h0, HarmonicConstantPointers& d0, const size_t N, const size_t batch_size){
    CopyToDevice(p,h,N,batch_size);
    void* kernelArgs[] = {(void*)&p, (void*)&d0};
    cudaLaunchCooperativeKernel((void*)kernel_HarmonicConstants<Block_Size_Pow_2>, dim3(batch_size,1,1), dim3(N,1,1), kernelArgs, sizeof(coord3d)*3*N, NULL);
    cudaDeviceSynchronize();
    cudaMemcpy(h0.bonds_0, d0.bonds_0, sizeof(coord3d)*N*batch_size , cudaMemcpyDeviceToHost);
    cudaMemcpy(h0.angles_0, d0.angles_0, sizeof(coord3d)*N*batch_size , cudaMemcpyDeviceToHost);
    cudaMemcpy(h0.outer_angles_m0, d0.outer_angles_m0, sizeof(coord3d)*N*batch_size , cudaMemcpyDeviceToHost);
    cudaMemcpy(h0.outer_angles_p0, d0.outer_angles_p0, sizeof(coord3d)*N*batch_size , cudaMemcpyDeviceToHost);
    cudaMemcpy(h0.angles_0, d0.angles_0, sizeof(coord3d)*N*batch_size , cudaMemcpyDeviceToHost);
    cudaMemcpy(h0.dihedrals_0, d0.dihedrals_0, sizeof(coord3d)*N*batch_size , cudaMemcpyDeviceToHost);
    cudaMemcpy(h0.outer_dihedrals_a0, d0.outer_dihedrals_a0, sizeof(coord3d)*N*batch_size , cudaMemcpyDeviceToHost);
    cudaMemcpy(h0.outer_dihedrals_m0, d0.outer_dihedrals_m0, sizeof(coord3d)*N*batch_size , cudaMemcpyDeviceToHost);
    cudaMemcpy(h0.outer_dihedrals_p0, d0.outer_dihedrals_p0, sizeof(coord3d)*N*batch_size , cudaMemcpyDeviceToHost);
}

template <int Block_Size_Pow_2>
void Gradients(DevicePointers& p, const HostPointers& h, const size_t N, const size_t batch_size, real_t* gradients){
    CopyToDevice(p,h,N,batch_size);
    void* kernelArgs[] = {(void*)&p};
    cudaLaunchCooperativeKernel((void*)kernel_Gradients<Block_Size_Pow_2>, dim3(batch_size,1,1), dim3(N,1,1), kernelArgs, sizeof(coord3d)*3*N, NULL);
    cudaDeviceSynchronize();
    cudaMemcpy(gradients, p.gradients, sizeof(coord3d)*N*batch_size , cudaMemcpyDeviceToHost);
}

template <int Block_Size_Pow_2>
void OptimizeBatch(DevicePointers& p, HostPointers& h, const size_t N, const size_t batch_size, const size_t MaxIter){
    getLastCudaError("One or more Mallocs Failed! \n");
    cudaMemcpy(p.X, h.h_X, sizeof(coord3d)*N*batch_size , cudaMemcpyHostToDevice);
    cudaMemcpy(p.neighbours, h.h_cubic_neighbours, sizeof(node_t)*3*N*batch_size, cudaMemcpyHostToDevice);
    cudaMemcpy(p.next_on_face, h.h_next_on_face, sizeof(node_t)*3*N*batch_size, cudaMemcpyHostToDevice);
    cudaMemcpy(p.prev_on_face, h.h_prev_on_face, sizeof(node_t)*3*N*batch_size, cudaMemcpyHostToDevice);
    cudaMemcpy(p.face_right, h.h_face_right, sizeof(uint8_t)*3*N*batch_size, cudaMemcpyHostToDevice);
    getLastCudaError("Memcpy Failed! \n");
    auto start = std::chrono::system_clock::now();
    void* kernelArgs[] = {
    (void*)&p,
    (void*)&N,
    (void*)&MaxIter
    };

    cudaLaunchCooperativeKernel((void*)kernel_OptimizeBatch<Block_Size_Pow_2>, dim3(batch_size, 1, 1), dim3(N, 1, 1), kernelArgs, sizeof(coord3d)*3*(Block_Size_Pow_2+1) + sizeof(real_t)*N, NULL);
    cudaDeviceSynchronize();
    auto end = std::chrono::system_clock::now();
    getLastCudaError("Failed to launch kernel: ");
    
    cudaMemcpy(h.h_X, p.X, sizeof(coord3d)*N*batch_size , cudaMemcpyDeviceToHost);
    getLastCudaError("Failed to copy back: ");
    
    std::cout << "Elapsed time: " << (end-start)/ 1ms << "ms\n" ;
    std::cout << "Estimated Performance " << ((real_t)(batch_size)/(std::chrono::duration_cast<std::chrono::microseconds>(end-start)).count()) * 1.0e6 << "Fullerenes/s \n";
}
};

