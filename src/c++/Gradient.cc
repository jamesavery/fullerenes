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
    const CubicArcs<N> neighbours;
    array<coord3d, N> X; // Current nucleus positions
    const array<uint8_t, N * 3> face_left, face_right;
    const array<node_t, N * 3> next_on_face, prev_on_face;

    // Optimal parameters. TODO: Fill out with the correct constants as constexpr's.
    constexpr array<real_t, 2> optimal_corner_cos_angles = {cos(M_PI * 108 / 180), cos(M_PI * 120 / 180)}; // indexed by face size: {p,h}  {0,1}
    constexpr array<real_t, 3> optimal_bond_lengths = {{1.479, 1.458}, {1.458, 1.401}}; // indexed by incident faces: {{pp,ph},{hp,pp}}
    constexpr array<real_t, 4> optimal_dih_cos_angles = {{{0.652358, 0.509674}, {0.509674, 0.345123}}, {{0.509674, 0.417884}, {0.417884, 0}}}; // indexed by ... faces: {{{ppp,pph},{php,phh}}, {{hpp,hph},{hhp,hhh}}}

    constexpr array<real_t, 2> angle_forces = {100.0, 100.0}; // {p,h}  {}
    constexpr array<real_t, 3> bond_forces = {{260.0, 390.0}, {390.0, 450.0}}; // {{pp,hp},{hp,hh}}
    constexpr array<real_t, 4> dih_forces{{{35.0, 65.0}, {65.0, 85.0}}, {{65.0, 85.0}, {85.0, 270.0}}}; // {{{ppp,pph},{pph,phh}},{{pph,phh},{phh,hhh}}}

    FullereneForcefieldEnergy(const CubicArcs<N> &neighbours, const array<coord3d, N> &X, const array<uint8_t, N * 3> &face_left, const array<uint8_t, N * 3> &face_right, const array<node_t, N * 3> &next_on_face, const array<node_t, N * 3> &prev_on_face) : 
        neighbours(neighbours), X(X), face_left(face_left), face_right(face_right), next_on_face(next_on_face), prev_on_face(prev_on_face) {}

    static inline coord3d bond_gradient(const real_t force_const, const real_t r0, const real_t rab, const coord3d ab_hat){
        return - force_const * (rab - r0) * ab_hat;
    }

    static inline coord3d angle_gradient(const real_t force_const, const real_t rab, const real_t rac, const real_t ang0, const coord3d &ab_hat, const coord3d &ac_hat)
    {
        real_t cos_angle = dot(ab_hat, ac_hat);
        coord3d grad = cos_angle * (ab_hat / rab + ac_hat / rac) - ab_hat / rac - ac_hat / rab;
        return force_const * (cos_angle - ang0) * grad;
    }

    static inline coord3d outer_angle_gradient(const real_t force_const, const real_t ang0, const real_t rab, coord3d &ab_hat, const coord3d &bp)
    {
        real_t rbp = bond_length(bp);
        coord3d bp_hat = bp / rbp;
        real_t cos_angle = dot(ab_hat, bp_hat);
        coord3d grad = (bp_hat - ab_hat * cos_angle) * rab;
        return force_const * (cos_angle - ang0) * grad;
    }

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

    coord3d gradient()
    {
        coord3d gradient = {0, 0, 0};

        for (node_t a = 0; a < N; a++)
        {
            array<node_t, 3> na = {neighbours[a * 3], neighbours[a * 3 + 1], neighbours[a * 3 + 2]}; // Neighbours b,c,d to a

            coord3d Xa = X[a];

            array<uint8_t, 3> f_r = {face_right[a * 3] - 5, face_right[a * 3 + 1] - 5, face_right[a * 3 + 2] - 5};
            array<uint8_t, 3> f_l = {face_left[a * 3] - 5, face_left[a * 3 + 1] - 5, face_left[a * 3 + 2]};

            array<array<uint8_t,3>, 3> f_l_bcd = {{face_left[na[0] * 3]]}, face_left[na[0] * 3 + 1], face_left[na[0] * 3 + 2]}, 
                                        {face_left[na[1] * 3], face_left[na[1] * 3 + 1], face_left[na[1] * 3 + 2]},
                                        {face_left[na[2] * 3], face_left[na[2] * 3 + 1], face_left[na[2] * 3 + 2]}};
            
            array<uint8_t,3> dihedral_face_sum_bcd = {sum(f_l_bcd[0]), sum(f_l_bcd[1]), sum(f_l_bcd[2])};

            array<coord3d, 3> Xbcd = {X[na[0]], X[na[1]], X[na[2]]};
            array<coord3d, 3> Xbp = {X[next_on_face[na[0]]], X[next_on_face[na[1]]], X[next_on_face[na[2]]]};
            array<coord3d, 3> Xbm = {X[prev_on_face[na[0]]], X[prev_on_face[na[1]]], X[prev_on_face[na[2]]]};

            array<real_t, 3> rab, rbp, rbm;

            array<coord3d, 3> abs, bps, bms, ab_hats, bp_hats, bm_hats;

            for (size_t i = 0; i < 3; i++)
            {
                tie(rab[i], abs[i]) = split_norm(Xbcd[i] - Xa);
                tie(rbp[i], bps[i]) = split_norm(Xbp[i] - Xa);
                tie(rbm[i], bms[i]) = split_norm(Xbm[i] - Xa);
                ab_hats[i] = abs[i] / rab[i];
                bp_hats[i] = bps[i] / rbp[i];
                bm_hats[i] = bms[i] / rbm[i];
            }
            for (size_t j = 0; j < 3; j++)
            {
                real_t r0 = optimal_bond_lengths[ f_l[j] + f_r[j] ];
                real_t ang0 = optimal_corner_cos_angles[ f_l[j] ];
                real_t dih0 = optimal_dih_cos_angles[ f_l[0] + f_l[1] + f_l[2]];
                real_t outer_ang0_m = optimal_corner_cos_angles[ f_r[j] ];
                real_t outer_ang0_p = optimal_corner_cos_angles[ f_l[j] ];
                real_t outer_dih0 = optimal_dih_cos_angles[ dihedral_face_sum_bcd[j] ];
                real_t bond_fc = bond_forces[ f_l[j] + f_r[j] ];
                real_t ang_fc = angle_forces[ f_l[j] ];
                real_t dih_fc = dih_forces[ f_l[0] + f_l[1] + f_l[2] ];
                real_t out_ang_fc_m = angle_forces[ f_r[j] ];
                real_t out_ang_fc_p = angle_forces[ f_l[j] ];
                real_t out_dih_fc = dih_forces[ dihedral_face_sum_bcd[j] ];

                coord3d bond_grad = bond_gradient(bond_fc, r0, rab[j], ab_hats[j]);
                coord3d angle_grad = angle_gradient(ang_fc, rab[j], rab[(j+1)%3], ang0, ab_hats[j], ab_hats[(j+1)%3]);
                coord3d dihedral_grad = dihedral_gradient(dih_fc, rab[j], dih0, ab_hats[j], abs[j], abs[(j+1)%3], abs[(j+2)%3]);

                coord3d outer_ang_grad_m = outer_angle_gradient(out_ang_fc_m, outer_ang0_m, rab[j], ab_hats[j], bms[j]);
                coord3d outer_ang_grad_p = outer_angle_gradient(out_ang_fc_p, outer_ang0_p, rab[j], ab_hats[j], bps[j]);

                coord3d outer_dihedral_a_grad = outer_dih_a_gradient(out_dih_fc, outer_dih0, rab[j], abs[j], bms[j], bps[j], ab_hats[j], bp_hats[j], bm_hats[j]);
                coord3d outer_dihedral_m_grad = outer_dih_m_gradient(out_dih_fc, outer_dih0, rab[j], abs[j], bms[j], bps[j], ab_hats[j], bp_hats[j], bm_hats[j]);
                coord3d outer_dihedral_p_grad = outer_dih_p_gradient(out_dih_fc, outer_dih0, rab[j], abs[j], bms[j], bps[j], ab_hats[j], bp_hats[j], bm_hats[j]);

                gradient =  gradient + bond_grad + angle_grad + dihedral_grad + outer_ang_grad_m + outer_ang_grad_p + outer_dihedral_a_grad + outer_dihedral_m_grad + outer_dihedral_p_grad;
                
            }
            

        }

        return gradient;
    }
};

int main()
{

    coord3d ab = {1, 1, 0.1};
    coord3d bm = {0, 1, 0.1};
    coord3d bp = {1, 0, 0.1};
    coord3d ab_hat = ab / bond_length(ab);
    coord3d bm_hat = bm / bond_length(bm);
    coord3d bp_hat = bp / bond_length(bp);

    coord3d result = FullereneForcefieldEnergy<60>::outer_dih_p_gradient(1, 2, bond_length(ab), ab, bm, bp, ab_hat, bp_hat, bm_hat);

    print_coord(result);
}