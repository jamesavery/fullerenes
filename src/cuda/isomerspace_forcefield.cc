
#include "coord3d.cc"
#include "forcefield_structs.cu"
#include <stdio.h>
#include "coord3d.cc"

struct Forcefield {
    
    size_t N;
    const NodeNeighbours* node_graph;
    const Constants* constants;

    struct ArcData
    {
        unsigned char j;
        size_t node;

        ArcData(const unsigned char node, const coord3d* X, const NodeNeighbours G){
            real_t r_rmp;
            coord3d ap, am, ab, ac, ad, mp;
            coord3d X_a = X[node]; coord3d X_b = X[G.neighbours[j]];
            //printf("Index: %d \n", a*3 + j);

            //Compute the arcs ab, ac, ad, bp, bm, ap, am, mp, bc and cd
            ab = (X_b - X_a);  r_rab = bond_length(ab); ab_hat = r_rab * ab;
            ac = (X[G.neighbours[j]%3] - X_a); r_rac = bond_length(ac); ac_hat = r_rac * ac; rab = non_resciprocal_bond_length(ab);
            ad = (X[G.neighbours[j]%3] - X_a); r_rad = bond_length(ad); ad_hat = r_rad * ad;
            
            coord3d bp = (X[G.next_on_face[j]] - X_b); bp_hat = unit_vector(bp);
            coord3d bm = (X[G.prev_on_face[j]] - X_b); bm_hat = unit_vector(bm);

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
    
    //3 FLOPs
    inline real_t harmonic_energy(const real_t p0, const real_t p) const{
        return (real_t)0.5*(p-p0)*(p-p0);
    }
    //4 FLOPs
    inline coord3d  harmonic_energy_gradient(const real_t p0, const real_t p, const coord3d gradp) const{
        return (p-p0)*gradp;     
    }
    //1 FLOP
    inline real_t bond() const {return rab;}
    //5 FLOPs
    inline real_t angle() const {return dot(ab_hat,ac_hat);}

    //Returns outer angle m, used only diagnostically.
    inline real_t outer_angle_m() const {return -dot(ab_hat, bm_hat);} //Compute outer angle. ab,bm

    //Returns outer angle p, used only diagnostically.
    inline real_t outer_angle_p() const{return -dot(ab_hat, bp_hat);} //Compute outer angle. ab,bp

    //Returns the inner dihedral angle for the current arc. Used here only for energy calculation, 
    //otherwise embedded in dihedral computation because the planes and angles that make up the dihedral angle computation are required for derivative computation.
    //50 FLOPs
    inline real_t dihedral() const 
    { 
        coord3d nabc, nbcd; real_t cos_b, cos_c, r_sin_b, r_sin_c;
        cos_b = dot(ba_hat,bc_hat);  r_sin_b = (real_t)1.0/sqrt((real_t)1.0 - cos_b*cos_b); nabc = cross(ba_hat, bc_hat) * r_sin_b;
        cos_c = dot(-bc_hat,cd_hat); r_sin_c = (real_t)1.0/sqrt((real_t)1.0 - cos_c*cos_c); nbcd = cross(-bc_hat,cd_hat) * r_sin_c;
        return dot(nabc, nbcd);
    }
    //Returns the Outer-dihedral-a wrt. current arc, only accessed diagnostically (internal coordinate).
    inline real_t outer_dihedral_a() const
    {
        coord3d nbam_hat, namp_hat; real_t cos_a, cos_m, r_sin_a, r_sin_m;
        cos_a = dot(ab_hat,am_hat); r_sin_a = (real_t)1.0/sqrt((real_t)1.0 - cos_a*cos_a); nbam_hat = cross(ab_hat,am_hat) * r_sin_a;
        cos_m = dot(-am_hat,mp_hat); r_sin_m = (real_t)1.0/sqrt((real_t)1.0 - cos_m*cos_m); namp_hat = cross(-am_hat,mp_hat) * r_sin_m;
        real_t cos_beta = dot(nbam_hat, namp_hat); //Outer Dihedral angle bam, amp
        return cos_beta;
    }
    //Returns the Outer-dihedral-m wrt. current arc, only accessed diagnostically (internal coordinate).
    inline real_t outer_dihedral_m() const
    {
        coord3d nbmp_hat, nmpa_hat; real_t cos_m, cos_p, r_sin_m, r_sin_p;
        cos_m = dot(mb_hat,mp_hat);  r_sin_m = (real_t)1.0/sqrt((real_t)1.0 - cos_m*cos_m); nbmp_hat = cross(mb_hat,mp_hat) * r_sin_m;
        cos_p = dot(-mp_hat,pa_hat); r_sin_p = (real_t)1.0/sqrt((real_t)1.0 - cos_p*cos_p); nmpa_hat = cross(-mp_hat,pa_hat) * r_sin_p;
        //Cosine to the outer dihedral angle constituted by the planes bmp and mpa
        real_t cos_beta = dot(nbmp_hat, nmpa_hat); //Outer dihedral angle bmp,mpa.
        return cos_beta;    
    }
    //Returns the Outer-dihedral-p wrt. current arc, only accessed diagnostically (internal coordinate).
    inline real_t outer_dihedral_p() const
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
    inline coord3d inner_angle_gradient(const Constants& c) const
    {   
        real_t cos_angle = angle(); //Inner angle of arcs ab,ac.
        coord3d grad = cos_angle * (ab_hat * r_rab + ac_hat * r_rac) - ab_hat * r_rac - ac_hat* r_rab; //Derivative of inner angle: Eq. 21. 
        return c.f_inner_angle[j] * harmonic_energy_gradient(c.angle0[j], cos_angle, grad); //Harmonic Energy Gradient: Eq. 21. multiplied by harmonic term.
    }
    //Computes gradient related to bending of outer angles. ~20 FLOPs
    inline coord3d outer_angle_gradient_m(const Constants& c) const
    {
        real_t cos_angle = -dot(ab_hat, bm_hat); //Compute outer angle. ab,bm
        coord3d grad = (bm_hat + ab_hat * cos_angle) * r_rab; //Derivative of outer angles Eq. 30. Buster Thesis
        return c.f_outer_angle_m[j] * harmonic_energy_gradient(c.outer_angle_m0[j],cos_angle,grad); //Harmonic Energy Gradient: Eq. 30 multiplied by harmonic term.
    }
    inline coord3d outer_angle_gradient_p(const Constants& c) const
    {
        real_t cos_angle = -dot(ab_hat, bp_hat); //Compute outer angle. ab,bp
        coord3d grad = (bp_hat + ab_hat * cos_angle) * r_rab; //Derivative of outer angles Eq. 28. Buster Thesis
        return c.f_outer_angle_p[j] * harmonic_energy_gradient(c.outer_angle_p0[j],cos_angle,grad); //Harmonic Energy Gradient: Eq. 28 multiplied by harmonic term.
    }
    // Chain rule terms for dihedral calculation
    //Computes gradient related to dihedral/out-of-plane term. ~75 FLOPs
    inline coord3d inner_dihedral_gradient(const Constants& c) const
    {
        coord3d nabc, nbcd; real_t cos_b, cos_c, r_sin_b, r_sin_c;
        cos_b = dot(ba_hat,bc_hat); r_sin_b = (real_t)1.0/sqrt((real_t)1.0 - cos_b*cos_b); nabc = cross(ba_hat, bc_hat) * r_sin_b;
        cos_c = dot(-bc_hat,cd_hat); r_sin_c = (real_t)1.0/sqrt((real_t)1.0 - cos_c*cos_c); nbcd = cross(-bc_hat,cd_hat) * r_sin_c;

        real_t cos_beta = dot(nabc, nbcd); //Inner dihedral angle from planes abc,bcd.
        real_t cot_b = cos_b * r_sin_b * r_sin_b; //cos(b)/sin(b)^2

        //Derivative w.r.t. inner dihedral angle F and G in Eq. 26
        coord3d grad = cross(bc_hat, nbcd) * r_sin_b * r_rab - ba_hat * cos_beta * r_rab + (cot_b * cos_beta * r_rab) * (bc_hat - ba_hat * cos_b);

        return c.f_inner_dihedral[j] * harmonic_energy_gradient(c.inner_dih0[j], cos_beta, grad); //Eq. 26.
    }

    //Computes gradient from dihedral angles constituted by the planes bam, amp ~162 FLOPs
    inline coord3d outer_dihedral_gradient_a(const Constants& c) const
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

        return c.f_outer_dihedral[j] * harmonic_energy_gradient(c.outer_dih0_a[j], cos_beta, grad);
    }

    //Computes gradient from dihedral angles constituted by the planes nbmp, nmpa ~92 FLOPs
    inline coord3d outer_dihedral_gradient_m(const Constants& c) const
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
        return c.f_outer_dihedral[j] * harmonic_energy_gradient(c.outer_dih0_m[j], cos_beta, grad);
    }

    //Computes gradient from dihedral angles constituted by the planes bpa, pam ~162 FLOPs
    inline coord3d outer_dihedral_gradient_p(const Constants& c) const
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
        return c.f_outer_dihedral[j] * harmonic_energy_gradient(c.outer_dih0_p[j], cos_beta, grad);
    }
    // Internal coordinate gradients
    inline coord3d bond_length_gradient(const Constants& c) const { return c.f_bond[j] * harmonic_energy_gradient(bond(),c.r0[j],ab_hat);}
    //Sum of angular gradient components.
    inline coord3d angle_gradient(const Constants& c) const { return inner_angle_gradient(c) + outer_angle_gradient_p(c) + outer_angle_gradient_m(c);}
    //Sum of inner and outer dihedral gradient components.
    inline coord3d dihedral_gradient(const Constants& c) const { return inner_dihedral_gradient(c) + outer_dihedral_gradient_a(c) + outer_dihedral_gradient_m(c) + outer_dihedral_gradient_p(c);}
    //coord3d flatness()             const { return ;  }   
    
    inline real_t bond_energy(const Constants& c) const {return (real_t)0.5 *c.f_bond[j] *harmonic_energy(bond(),c.r0[j]);}
    inline real_t bend_energy(const Constants& c) const {return c.f_inner_angle[j]* harmonic_energy(angle(),c.angle0[j]);}
    inline real_t dihedral_energy(const Constants& c) const {return c.f_inner_dihedral[j]* harmonic_energy(dihedral(),c.inner_dih0[j]);}
    //Harmonic energy contribution from bond stretching, angular bending and dihedral angle bending.
    //71 FLOPs
    inline real_t energy(const Constants& c) const {return bond_energy(c) + bend_energy(c) + dihedral_energy(c); }
    //Sum of bond, angular and dihedral gradient components.
    inline coord3d gradient(const Constants& c) const{return bond_length_gradient(c) + angle_gradient(c) + dihedral_gradient(c);}

    inline real_t GSS(coord3d* X, coord3d& r0, coord3d* X1, coord3d* X2) const{
        constexpr real_t tau = (real_t)0.6180339887;
        //Line search x - values;
        real_t a = 0.0; real_t b = (real_t)1.0;
        
        real_t x1,  x2;
        x1 = (a + (1 - tau) * (b - a));
        x2 = (a + tau * (b - a));
        //Actual coordinates resulting from each traversal 
        X1[node] = X[node] + x1 * r0;
        X2[node] = X[node] + x2 * r0;
        real_t f1 = energy(X1);
        real_t f2 = energy(X2);

        for (uint8_t i = 0; i < 50; i++){
            if (f1 > f2){
                a = x1;
                x1 = x2;
                f1 = f2;
                x2 = a + tau * (b - a);
                X2[node] = X[node] + x2 * r0;
                f2 = energy(X2);
            }else
            {
                b = x2;
                x2 = x1;
                f2 = f1;
                x1 = a + ((real_t)1.0 - tau) * (b - a);
                X1[node] = X[node] + x1 * r0;
                f1 = energy(X1);
            }
        }
        if (f1 > energy(X)) {return (real_t)0.0;}
        //Line search coefficient
        real_t alpha = (a+b)/(real_t)2.0;
        return alpha;
    }
    };

    inline  void CG(coord3d* X, coord3d* X1, coord3d* X2, const size_t MaxIter)
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
            alpha = LINESEARCH_METHOD(X,s,X1,X2);
            if (alpha > (real_t)0.0){X1[node] = X[node] + alpha * s;}
            g1 = gradient(X1);
            //Polak Ribiere method
            g0_norm2 = reduction(sdata, dot(g0, g0));
            beta = max(reduction(sdata, dot(g1, (g1 - g0))) / g0_norm2,(real_t)0.0);

            if (alpha > (real_t)0.0){X[node] = X1[node];}else{ g1 = g0; beta = (real_t) 0.0;}
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