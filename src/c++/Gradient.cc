#include <chrono>
#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
using namespace std;

#include <array>
#include <tuple>

#include "coord3d.cc"

typedef array<real_t, 3> coord3d;
typedef uint16_t node_t; // NB: Increase to int32 to do more than 65k atoms

inline void print_real(real_t a)
{
    cout << a << "\n";
}

inline uint8_t sum(array<uint8_t,3> a){
    return a[0] + a[1] + a[2];
}

template <int N>
using CubicArcs = array<node_t, N * 3>;

template <int N>
class FullereneForcefieldEnergy
{
public:
    const CubicArcs<N> neighbours; //Bookkeeping array for storing the indices of neighbour nodes: b,c,d
    array<coord3d, N> X; // Current nucleus positions

    const array<uint8_t, N * 3> face_right; //Faces to the right of the edges ab, ac, ad
    const array<node_t, N * 3> next_on_face, prev_on_face; //Next node on the face to the left of the edges ab, ac, ad

    
    //All parameter arrays are indexed by a binary sum, 0,1,2,3,4,...
    //Pentagons = 0
    //Hexagons = 1
    //PPP = 0, {HPP, PHP, PPH} = 1, {PHH, HPH, HHP} = 2, {HHH} = 3
    const array<real_t, 2> optimal_corner_cos_angles = {cos(M_PI * 108 / 180), cos(M_PI * 120 / 180)}; 
    const array<real_t, 3> optimal_bond_lengths = {1.479, 1.458, 1.401}; 
    const array<real_t, 4> optimal_dih_cos_angles = {cos(0.652358), cos(0.509674), cos(0.417884), cos(0)}; 

    const array<real_t, 2> angle_forces = {100.0, 100.0}; 
    const array<real_t, 3> bond_forces = {260.0, 390.0, 450.0}; 
    const array<real_t, 4> dih_forces = {35.0, 65.0, 85.0, 270.0}; 

    FullereneForcefieldEnergy(const CubicArcs<N> &neighbours, const array<coord3d, N> &X, const array<uint8_t, N * 3> &face_right, const array<node_t, N * 3> &next_on_face, const array<node_t, N * 3> &prev_on_face) : 
        neighbours(neighbours), X(X), face_right(face_right), next_on_face(next_on_face), prev_on_face(prev_on_face) {}

    //Computes gradient related to stretch term.
    static inline coord3d bond_gradient(const real_t force_const, const real_t r0, const real_t rab, const coord3d ab_hat){
        return - force_const * (rab - r0) * ab_hat;
    }

    //Computes gradient related to bending term.
    static inline coord3d angle_gradient(const real_t force_const, const real_t rab, const real_t rac, const real_t ang0, const coord3d &ab_hat, const coord3d &ac_hat)
    {   
        real_t cos_angle = dot(ab_hat, ac_hat);
        coord3d grad = cos_angle * (ab_hat / rab + ac_hat / rac) - ab_hat / rac - ac_hat / rab;
        return force_const * (cos_angle - ang0) * grad;
    }

    //Computes gradient related to bending of angles.
    static inline coord3d outer_angle_gradient(const real_t force_const, const real_t ang0, const real_t rab, coord3d &ab_hat, const coord3d &bp)
    {   
        coord3d ba_hat = -ab_hat;
        real_t rbp = bond_length(bp);
        coord3d bp_hat = bp / rbp;
        real_t cos_angle = dot(ba_hat, bp_hat);
        coord3d grad = (bp_hat - ba_hat * cos_angle) / rab;

        return force_const * (cos_angle - ang0) * grad;
    }

    //Computes gradient related to dihedral/out-of-plane term.
    static inline coord3d dihedral_gradient(const real_t force_const, const real_t rab, const real_t dih0, const coord3d &ab_hat, const coord3d &ab, const coord3d &ac, const coord3d &ad)
    {
        coord3d bc = ac - ab;
        coord3d cd = ad - ac;

        real_t rbc = bond_length(bc);
        real_t rcd = bond_length(cd);

        //Compute normalized vectors.
        coord3d bc_hat = bc / rbc;
        coord3d cd_hat = cd / rcd;
        coord3d cb_hat = -bc_hat;
        coord3d ba_hat = -ab_hat;

        real_t cos_b = dot(ba_hat, bc_hat);
        real_t cos_c = dot(cb_hat, cd_hat);
        real_t sin_b = sqrt(1 - cos_b * cos_b);
        real_t sin_c = sqrt(1 - cos_c * cos_c);

        coord3d nabc = cross(ba_hat, bc_hat) / sin_b;
        coord3d nbcd = cross(cb_hat, cd_hat) / sin_c;

        real_t cos_beta = dot(nabc, nbcd);
        real_t cot_b = cos_b / sin_b;

        coord3d d_dih_a = cross(bc_hat, nbcd) / (sin_b * rab) - (ba_hat * cos_beta) / rab + (cot_b * cos_beta) / (sin_b * rab) * (bc_hat - ba_hat * cos_b);

        return force_const * (cos_beta - dih0) * d_dih_a;
    }

    //Computes gradient from dihedral angles constituted by the planes nbam, namp
    static inline coord3d outer_dih_a_gradient(const real_t force_const, const real_t dih0, const real_t rab, const coord3d &ab, const coord3d &bm, const coord3d &bp, coord3d &ab_hat, coord3d &bp_hat, coord3d &bm_hat)
    {
        coord3d ba_hat = -ab_hat;
        coord3d ba = -ab;

        coord3d am = bm - ba;
        real_t ram = bond_length(am);
        coord3d am_hat = am / ram;
        coord3d ma_hat = -am_hat;
        coord3d mp = bp - bm;

        real_t rmp = bond_length(mp);

        coord3d mp_hat = mp / rmp;

        real_t cos_a = dot(ab_hat, am_hat);
        real_t cos_m = dot(ma_hat, mp_hat);

        real_t sin_a = sqrt(1 - cos_a * cos_a);
        real_t sin_m = sqrt(1 - cos_m * cos_m);

        coord3d nbam_hat = cross(ab_hat, am_hat) / sin_a;
        coord3d namp_hat = cross(ma_hat, mp_hat) / sin_m;

        real_t cos_beta = dot(nbam_hat, namp_hat);
        real_t cot_a = cos_a / sin_a;
        real_t cot_m = cos_m / sin_m;

        coord3d d_dih_a = (ab_hat * cos_beta / rab) - cross(am_hat, namp_hat) / (rab * sin_a) +
                          am_hat * cos_beta / ram - cross(namp_hat, ab_hat) / (ram * sin_a) +
                          (cot_a * cos_beta / sin_a) * (ab_hat * cos_a / rab - am_hat / rab + am_hat * cos_a / ram - ab_hat / ram) +
                          cross(mp_hat, nbam_hat) / (ram * sin_m) - ma_hat * cos_beta / ram +
                          (cot_m * cos_beta / sin_m) * (mp_hat / ram - ma_hat * cos_m / ram);

        return force_const * (cos_beta - dih0) * d_dih_a;
    }
    
    //Computes gradient from dihedral angles constituted by the planes nbmp, nmpa
    static inline coord3d outer_dih_m_gradient(const real_t force_const, const real_t dih0, const real_t rab, const coord3d &ab, const coord3d &bm, const coord3d &bp, coord3d &ab_hat, coord3d &bp_hat, coord3d &bm_hat)
    {
        coord3d ba_hat = -ab_hat;
        coord3d ba = -ab;

        coord3d mp = bp - bm;
        coord3d pa = ba - bp;

        real_t rmp = bond_length(mp);
        real_t rpa = bond_length(pa);

        coord3d mp_hat = mp / rmp;
        coord3d pa_hat = pa / rpa;
        coord3d mb_hat = -bm_hat;
        coord3d pm_hat = -mp_hat;

        real_t cos_m = dot(mb_hat, mp_hat);
        real_t cos_p = dot(pm_hat, pa_hat);
        real_t sin_m = sqrt(1 - cos_m * cos_m);
        real_t sin_p = sqrt(1 - cos_p * cos_p);

        coord3d nbmp_hat = cross(mb_hat, mp_hat) / sin_m;
        coord3d nmpa_hat = cross(pm_hat, pa_hat) / sin_p;

        real_t cos_beta = dot(nbmp_hat, nmpa_hat);
        real_t cot_p = cos_p / sin_p;

        coord3d d_dih_m = cross(nbmp_hat, pm_hat) / (rpa * sin_p) - pa_hat * cos_beta / rpa + (cot_p * cos_beta / (sin_p * rpa)) * (pm_hat - pa_hat * cos_p);

        return force_const * (cos_beta - dih0) * d_dih_m;
    }
    
    //Computes gradient from dihedral angles constituted by the planes nbpa, npam
    static inline coord3d outer_dih_p_gradient(const real_t force_const, const real_t dih0, const real_t rab, const coord3d &ab, const coord3d &bm, const coord3d &bp, coord3d &ab_hat, coord3d &bp_hat, coord3d &bm_hat)
    {
        coord3d pb_hat = -bp_hat;
        coord3d ba = -ab;
        coord3d pa = ba - bp;

        real_t rpa = bond_length(pa);

        coord3d pa_hat = pa / rpa;
        coord3d ap_hat = -pa_hat;
        coord3d am = bm - ba;

        real_t ram = bond_length(am);

        coord3d am_hat = am / ram;

        real_t cos_p = dot(pb_hat, pa_hat);
        real_t cos_a = dot(ap_hat, am_hat);

        real_t sin_p = sqrt(1 - cos_p * cos_p);
        real_t sin_a = sqrt(1 - cos_a * cos_a);

        coord3d nbpa_hat = cross(pb_hat, pa_hat) / sin_p;
        coord3d npam_hat = cross(ap_hat, am_hat) / sin_a;

        real_t cos_beta = dot(nbpa_hat, npam_hat);
        real_t cot_p = cos_p / sin_p;
        real_t cot_a = cos_a / sin_a;

        coord3d d_dih_p = cross(npam_hat, pb_hat) / (rpa * sin_p) - pa_hat * cos_beta / rpa +
                          (cot_p * cos_beta / (rpa * sin_p)) * (pb_hat - pa_hat * cos_p) +
                          ap_hat * cos_beta / rpa - cross(am_hat, nbpa_hat) / (rpa * sin_a) + am_hat * cos_beta / ram -
                          cross(nbpa_hat, ap_hat) / (ram * sin_a) + (cot_a * cos_beta / sin_a) * (ap_hat * cos_a / rpa - am_hat / rpa + am_hat * cos_a / ram - ap_hat / ram);
        return force_const * (cos_beta - dih0) * d_dih_p;
    }

    array<coord3d,N> gradient()
    {
        array<coord3d,N> gradient;
        for (node_t a = 0; a < N; a++) 
        {   
            //Fetching neighbour indices.
            array<node_t, 3> na = {neighbours[a * 3], neighbours[a * 3 + 1], neighbours[a * 3 + 2]}; // Neighbours b,c,d to a

            //Fetching current node coordinate.
            coord3d Xa = X[a];

            //Fetching face information to the left and right of the edges: ab, ac, ad.
            array<uint8_t, 3> f_r = {face_right[a * 3] - 5, face_right[a * 3 + 1] - 5, face_right[a * 3 + 2] - 5};
            array<uint8_t, 3> f_l = {face_right[a * 3 + 2] - 5, face_right[a * 3] - 5, face_right[a * 3 + 1] - 5};

            //Fetching face information to the right of the edges {ba, bm, bp}, {ca, cm, cp}, {da, dm, dp}
            array<array<uint8_t,3>, 3> f_r_bcd = {face_right[na[0] * 3]-5, face_right[na[0] * 3 + 1]-5, face_right[na[0] * 3 + 2]-5, 
                                        face_right[na[1] * 3]-5, face_right[na[1] * 3 + 1]-5, face_right[na[1] * 3 + 2]-5,
                                        face_right[na[2] * 3]-5, face_right[na[2] * 3 + 1]-5, face_right[na[2] * 3 + 2]-5};
            
            //Pre calculate the sum of the face information for each neighbour. Used to index dihedral constants.
            array<uint8_t,3> dihedral_face_sum_bcd = {sum(f_r_bcd[0]), sum(f_r_bcd[1]), sum(f_r_bcd[2])};

            //Fetching coordinates for neighbour nodes b,c,d.
            array<coord3d, 3> Xbcd = {X[na[0]], X[na[1]], X[na[2]]};

            //Fetching coordinates for outer nodes.
            array<coord3d, 3> Xbp = {X[next_on_face[a*3]], X[next_on_face[a*3 + 1]], X[next_on_face[a*3 + 2]]};
            array<coord3d, 3> Xbm = {X[prev_on_face[a*3]], X[prev_on_face[a*3 + 1]], X[prev_on_face[a*3 + 2]]};

            //Vector norms.
            array<real_t, 3> rab, rbp, rbm;

            //Common vectors.
            array<coord3d, 3> abs, bps, bms, ab_hats, bp_hats, bm_hats;

            coord3d node_gradient = {0,0,0};

            for (size_t i = 0; i < 3; i++)
            {
                tie(rab[i], abs[i]) = split_norm(Xbcd[i] - Xa);
                tie(rbp[i], bps[i]) = split_norm(Xbp[i] - Xbcd[i]);
                tie(rbm[i], bms[i]) = split_norm(Xbm[i] - Xbcd[i]);

                ab_hats[i] = abs[i] / rab[i];
                bp_hats[i] = bps[i] / rbp[i];
                bm_hats[i] = bms[i] / rbm[i];
            }
            for (size_t j = 0; j < 3; j++)
            {   
                //Indexing of 'Equillibrium' parameters and related force constants, through face information.
                real_t r0 = optimal_bond_lengths[ f_l[j] + f_r[j] ];
                real_t ang0 = optimal_corner_cos_angles[ f_r[j] ];
                real_t dih0 = optimal_dih_cos_angles[ f_l[0] + f_l[1] + f_l[2]];
                real_t outer_ang0_m = optimal_corner_cos_angles[ f_l[j] ];
                real_t outer_ang0_p = optimal_corner_cos_angles[ f_r[j] ];
                real_t outer_dih0 = optimal_dih_cos_angles[ dihedral_face_sum_bcd[j] ];
                real_t bond_fc = bond_forces[ f_l[j] + f_r[j] ];
                real_t ang_fc = angle_forces[ f_l[j] ];
                real_t dih_fc = dih_forces[ f_l[0] + f_l[1] + f_l[2] ];
                real_t out_ang_fc_m = angle_forces[ f_r[j] ];
                real_t out_ang_fc_p = angle_forces[ f_l[j] ];
                real_t out_dih_fc = dih_forces[ dihedral_face_sum_bcd[j] ];

                //Computation of gradient terms.
                coord3d bond_grad = bond_gradient(bond_fc, r0, rab[j], ab_hats[j]);
                coord3d angle_grad = angle_gradient(ang_fc, rab[j], rab[(j+1)%3], ang0, ab_hats[j], ab_hats[(j+1)%3]);
                coord3d dihedral_grad = dihedral_gradient(dih_fc, rab[j], dih0, ab_hats[j], abs[j], abs[(j+1)%3], abs[(j+2)%3]);

                coord3d outer_ang_grad_m = outer_angle_gradient(out_ang_fc_m, outer_ang0_m, rab[j], ab_hats[j], bms[j]);
                coord3d outer_ang_grad_p = outer_angle_gradient(out_ang_fc_p, outer_ang0_p, rab[j], ab_hats[j], bps[j]);
                
                coord3d outer_dihedral_a_grad = outer_dih_a_gradient(out_dih_fc, outer_dih0, rab[j], abs[j], bms[j], bps[j], ab_hats[j], bp_hats[j], bm_hats[j]);
                coord3d outer_dihedral_m_grad = outer_dih_m_gradient(out_dih_fc, outer_dih0, rab[j], abs[j], bms[j], bps[j], ab_hats[j], bp_hats[j], bm_hats[j]);
                coord3d outer_dihedral_p_grad = outer_dih_p_gradient(out_dih_fc, outer_dih0, rab[j], abs[j], bms[j], bps[j], ab_hats[j], bp_hats[j], bm_hats[j]);

                node_gradient = node_gradient + bond_grad + angle_grad + dihedral_grad + outer_ang_grad_m + outer_ang_grad_p + outer_dihedral_a_grad + outer_dihedral_m_grad + outer_dihedral_p_grad; 

            }
            gradient[a] = gradient[a] + node_gradient;
        }
        return gradient;
    }
};

int main()
{
    const int size = 60;

    //Data necessary to compute gradients for the C60ih Fullerene
    CubicArcs<size> cubic_neighbours  = {4,15,1,0,12,2,1,9,3,2,5,4,3,8,0,3,11,6,5,22,7,6,20,8,7,18,4,2,14,10,9,26,11,10,25,5,1,17,13,12,30,14,13,29,9,0,19,16,15,34,17,16,33,12,8,21,19,18,37,15,7,24,21,20,38,18,6,25,23,22,42,24,23,40,20,11,28,22,10,29,27,26,45,28,27,44,25,14,32,26,13,33,31,30,48,32,31,47,29,17,36,30,16,37,35,34,51,36,35,50,33,19,39,34,21,41,39,38,53,37,24,43,41,40,54,38,23,44,43,42,55,40,28,46,42,27,47,46,45,57,44,32,49,45,31,50,49,48,58,47,36,52,48,35,53,52,51,59,50,39,54,51,41,56,53,43,57,56,55,59,54,46,58,55,49,59,57,52,56,58};    
    array<node_t, size*3> next_on_face = {8,16,2,15,13,3,12,10,4,9,6,0,5,18,1,2,25,7,11,23,8,22,21,4,20,19,3,1,29,11,14,27,5,26,22,3,0,33,14,17,31,9,30,26,2,4,37,17,19,35,12,34,30,1,7,38,15,21,34,0,6,40,18,24,39,8,5,28,24,25,43,20,42,41,7,10,44,6,9,32,28,29,46,25,45,42,11,13,47,10,12,36,32,33,49,29,48,45,14,16,50,13,15,39,36,37,52,33,51,48,17,18,53,16,20,54,37,41,51,19,23,55,38,43,53,21,22,46,40,44,56,24,27,57,23,26,49,44,47,55,28,31,58,27,30,52,47,50,57,32,35,59,31,34,54,50,53,58,36,38,56,35,40,59,39,42,58,54,57,52,41,45,59,43,48,56,46,51,55,49};
    array<node_t, size*3> prev_on_face = {3,19,12,4,17,9,0,14,5,1,11,8,2,7,15,4,10,22,3,25,20,5,24,18,6,21,0,3,13,26,2,29,25,9,28,6,2,16,30,1,33,29,12,32,10,1,18,34,0,37,33,15,36,13,4,20,37,8,39,16,8,23,38,7,41,19,7,11,42,6,44,40,22,43,21,5,27,23,11,14,45,10,47,44,26,46,22,9,31,27,14,17,48,13,50,47,30,49,26,12,35,31,17,19,51,16,53,50,34,52,30,15,38,35,18,40,53,21,54,34,20,42,54,24,56,39,24,28,55,23,57,41,25,45,43,28,32,57,27,58,42,29,48,46,32,36,58,31,59,45,33,51,49,36,39,59,35,56,48,37,41,52,38,55,51,40,46,59,43,58,53,44,49,56,47,52,55,50,54,57};
    array<uint8_t, size*3> face_right   = {6,6,5,6,6,5,6,6,5,6,6,5,6,6,5,6,5,6,5,6,6,6,5,6,5,6,6,6,5,6,5,6,6,6,5,6,6,5,6,5,6,6,6,5,6,6,5,6,5,6,6,6,5,6,5,6,6,6,5,6,6,6,5,6,6,5,5,6,6,6,5,6,5,6,6,6,6,5,5,6,6,6,5,6,5,6,6,6,6,5,5,6,6,6,5,6,5,6,6,6,6,5,5,6,6,6,5,6,5,6,6,6,6,5,6,5,6,5,6,6,5,6,6,6,5,6,6,6,5,6,6,5,5,6,6,6,6,5,6,6,5,5,6,6,6,6,5,6,6,5,5,6,6,6,6,5,6,6,5,5,6,6,6,6,5,6,5,6,5,6,6,6,5,6,6,5,6,6,5,6};
    array<coord3d, size> X = {2.82489,8.43261e-17,14.2017,8.0219,0.604118,12.0396,10.4266,6.01981,8.04461,4.53413,6.64511,12.0396,1.41245,2.44643,14.2017,3.29447,11.5801,8.04461,-3.29447,11.5801,8.04461,-4.53413,6.64511,12.0396,-1.41245,2.44643,14.2017,12.299,7.10085,2.82489,9.58483,10.4795,-2.82489,4.63059,13.4256,2.82489,11.6759,-2.93695,8.04461,13.9422,-2.70257,2.82489,13.8679,3.06098,-2.82489,1.41245,-2.44643,14.2017,3.48777,-7.24923,12.0396,8.38143,-8.64315,8.04461,-2.82489,5.87089e-16,14.2017,-1.41245,-2.44643,14.2017,-10.4266,6.01981,8.04461,-8.0219,0.604118,12.0396,-4.63059,13.4256,2.82489,-9.58483,10.4795,-2.82489,-12.299,7.10085,2.82489,2.39962e-16,14.2017,-2.82489,8.97979,8.0197,-8.04461,3.58009,7.20408,-12.0396,1.07573e-16,12.0396,-8.04461,11.4352,3.76688,-8.04461,12.299,-7.10085,-2.82489,10.4266,-6.01981,-8.04461,8.02896,-0.501593,-12.0396,9.31158,-10.723,2.82489,-2.56576e-15,-12.0396,8.04461,-2.91345e-15,-14.2017,2.82489,4.28306,-13.5404,-2.82489,-3.48777,-7.24923,12.0396,-11.6759,-2.93695,8.04461,-8.38143,-8.64315,8.04461,-13.8679,3.06098,-2.82489,-13.9422,-2.70257,2.82489,-8.97979,8.0197,-8.04461,-11.4352,3.76688,-8.04461,-3.58009,7.20408,-12.0396,1.33536,2.48934,-14.2017,-1.33536,2.48934,-14.2017,2.82351,-0.088213,-14.2017,4.44887,-6.70249,-12.0396,1.48815,-2.40113,-14.2017,2.45537,-11.7866,-8.04461,-4.28306,-13.5404,-2.82489,-2.45537,-11.7866,-8.04461,-9.31158,-10.723,2.82489,-12.299,-7.10085,-2.82489,-8.02896,-0.501593,-12.0396,-10.4266,-6.01981,-8.04461,-2.82351,-0.088213,-14.2017,-1.48815,-2.40113,-14.2017,-4.44887,-6.70249,-12.0396};

    //Test gradient computation
    FullereneForcefieldEnergy<60> forcefield = FullereneForcefieldEnergy<60>(cubic_neighbours, X, face_right, next_on_face, prev_on_face);
    array<coord3d, size> result = forcefield.gradient();
    for (size_t i = 0; i < size; i++)
    {
        print_coord(result[i]);
    }
    
    
}