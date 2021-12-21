
#include "coord3d.cc"
#include "forcefield_structs.cu"
#include <stdio.h>
#include "coord3d.cc"

struct Forcefield {
    
    size_t N;
    const NodeGraph* node_graph;
    const Constants* constants;

    struct ArcData
    {
        unsigned char j;
        size_t node;

        ArcData(const unsigned char node, const coord3d* X, const NodeGraph G){
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
    };
    

}