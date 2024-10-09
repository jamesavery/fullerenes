#include <CL/sycl.hpp>
#include "numeric"
#include <vector>
#include <tuple>
#include <iterator>
#include <type_traits>
#include "fstream"
#define USE_MAX_NORM 0
#define SQRT cl::sycl::sqrt
#include <fullerenes/sycl-headers/hessian-kernel.hh>
#include "queue-impl.cc"
#include "forcefield-includes.cc"

template <ForcefieldType FFT, typename T, typename K>
struct ForceField
{
    TEMPLATE_TYPEDEFS(T, K);
    typedef Constants<T, K> Constants;
    typedef mat3<T> mat3;

    const NodeNeighbours<K> node_graph; // Contains face-information and neighbour-information. Both of which are constant in the lifespan of this struct.
    const Constants constants;          // Contains force-constants and equillibrium-parameters. Constant in the lifespan of this struct.

    size_t node_id;
    size_t N;
    const sycl::group<1> cta;
    real_t *sdata; // Pointer to start of L1 cache array, used exclusively for reduction.

    ForceField(const NodeNeighbours<K> &G,
               const Constants &c,
               sycl::group<1> cta,
               real_t *sdata) : node_graph(G), constants(c), cta(cta), sdata(sdata)
    {
        node_id = cta.get_local_linear_id();
        N = cta.get_local_linear_range();
    }

    struct FaceData
    {
        coord3d Xa;
        symMat3<T> A;
        coord3d n_f;     // normalized normal vector to face-plane
        real_t lambda_f; // Smallest eigenvalue defining the flatness of the face
        coord3d lambdas;
        coord3d centroid;
        node3 face_neighbours;
        node_t N;
        node_t node_id;
        sycl::group<1> cta;
        // 84 + 107 FLOPS
        FaceData(const sycl::group<1> &cta, const Span<coord3d> X, const NodeNeighbours<K> &G) : cta(cta)
        {
            node_id = cta.get_local_linear_id();
            N = cta.get_local_linear_range();

            face_neighbours = G.face_neighbours;
            Xa = X[node_id];
            // There are only N/2 + 2 faces. (Nf  =  N/2 + 1)
            if (node_id < N / 2 + 2)
            {
                coord3d Xf[6] = {X[G.face_nodes[0]], X[G.face_nodes[1]], X[G.face_nodes[2]], X[G.face_nodes[3]], X[G.face_nodes[4]]};
                // If pentagon set to 0 otherwise get the 6th node coordinates.
                if (G.face_size == 6)
                {
                    Xf[5] = X[G.face_nodes[5]];
                }
                else
                {
                    Xf[5] = {(real_t)0., (real_t)0., (real_t)0.};
                }
                centroid = (Xf[0] + Xf[1] + Xf[2] + Xf[3] + Xf[4] + Xf[5]) / (T)G.face_size;
                // Centralise coordinate system to centroid of the face
                Xf[0] -= centroid;
                Xf[1] -= centroid;
                Xf[2] -= centroid;
                Xf[3] -= centroid;
                Xf[4] -= centroid;
                if (G.face_size == 6)
                {
                    Xf[5] -= centroid;
                }
                auto a = Xf[0][0] * Xf[0][0] + Xf[1][0] * Xf[1][0] + Xf[2][0] * Xf[2][0] + Xf[3][0] * Xf[3][0] + Xf[4][0] * Xf[4][0] + Xf[5][0] * Xf[5][0],
                     b = Xf[0][0] * Xf[0][1] + Xf[1][0] * Xf[1][1] + Xf[2][0] * Xf[2][1] + Xf[3][0] * Xf[3][1] + Xf[4][0] * Xf[4][1] + Xf[5][0] * Xf[5][1],
                     c = Xf[0][0] * Xf[0][2] + Xf[1][0] * Xf[1][2] + Xf[2][0] * Xf[2][2] + Xf[3][0] * Xf[3][2] + Xf[4][0] * Xf[4][2] + Xf[5][0] * Xf[5][2],
                     d = Xf[0][1] * Xf[0][1] + Xf[1][1] * Xf[1][1] + Xf[2][1] * Xf[2][1] + Xf[3][1] * Xf[3][1] + Xf[4][1] * Xf[4][1] + Xf[5][1] * Xf[5][1],
                     e = Xf[0][1] * Xf[0][2] + Xf[1][1] * Xf[1][2] + Xf[2][1] * Xf[2][2] + Xf[3][1] * Xf[3][2] + Xf[4][1] * Xf[4][2] + Xf[5][1] * Xf[5][2],
                     f = Xf[0][2] * Xf[0][2] + Xf[1][2] * Xf[1][2] + Xf[2][2] * Xf[2][2] + Xf[3][2] * Xf[3][2] + Xf[4][2] * Xf[4][2] + Xf[5][2] * Xf[5][2];
                // Xf * Xf^FFT In closed form.
                A = symMat3(a, b, c, d, e, f);

                // A is positive-semi-definite so all eigenvalues are non-negative
                lambdas = A.eigenvalues();
                lambda_f = d_min(d_min(lambdas[0], lambdas[1]), lambdas[2]);
            }
        }
        // 3 FLOPs
        /**
         * Computes the harmonic energy contribution of one term.
         *
         * @param[in] p0 Equillibrium parameter
         * @param[in] p Current parameter
         * @return Hooke's law harmonic energy contribution of the term
         */
        real_t harmonic_energy(const real_t p0, const real_t p) const
        {
            return (real_t)0.5 * (p - p0) * (p - p0);
        }

        /** @brief Compute the flatness of the threadIdx^th face in the isomer
         *  @return The flatness of the threadIdx^th face in the isomer
         */
        real_t flatness() const { return node_id < N / 2 + 2 ? lambda_f : (real_t)0.; }

        // 4 FLOPs
        /**
         * Computes the harmonic energy gradient contribution of one term.
         *
         * @param[in] p0 Equillibrium parameter
         * @param[in] p Current parameter
         * @param[in] gradp Gradient of the parameter w.r.t. the particle position
         * @return Hooke's law harmonic energy gradient contribution of the term
         */
        coord3d harmonic_energy_gradient(const real_t p0, const real_t p, const coord3d gradp) const
        {
            return (p - p0) * gradp;
        }

        /**
         * @brief Compute the flatness energy contribution of the threadIdx^th face in the isomer.
         *
         * @param c The forcefield constants for the threadIdx^th node.
         * @return The flatness energy.
         */
        real_t flatness_energy(const Constants &c) const
        {
            return c.f_flat * harmonic_energy(flatness(), (real_t)0.);
        }

        /**
         * @brief Compute the gradient of the flatness w.r.t to the threadIdx^th atom in the isomer.
         * @param c The forcefield constants for the threadIdx^th node.
         * @param cache A pointer to a cache of minimum size Nf * 2 * sizeof(coord3d) bytes.
         * @return The flatness energy gradient.
         */
        coord3d flatness_gradient(const Constants &c, coord3d *cache) const
        {
            coord3d *centroids = reinterpret_cast<coord3d *>(cache);
            coord3d *norms = reinterpret_cast<coord3d *>(cache + N / 2 + 2);
            if (node_id < N / 2 + 2)
            {
                auto lam = lambda_f;
                centroids[node_id] = centroid;
                norms[node_id] = A.eigenvector3x3(lam);
            }
            sycl::group_barrier(cta);

            coord3d grad = {(real_t)0., (real_t)0., (real_t)0.};
            for (unsigned char j = 0; j < 3; j++)
                grad += dot(Xa - centroids[d_get(face_neighbours, j)], norms[d_get(face_neighbours, j)]) * norms[d_get(face_neighbours, j)];
            return c.f_flat * (real_t)2. * grad;
        }
    };

    // Container for all energy and gradient evaluations with respect to an arc, eg. AB, AC or AD.
    struct ArcData
    {
        // 124 FLOPs;
        int j;
        /**
         * @brief Construct a new ArcData object
         * @param j The index of the arc, eg. 0 for ab, 1 for ac and 2 for ad.
         * @param X The coordinates of all nodes in the isomer.
         * @param G The neighbour information for the threadIdx^th node.
         * @return A new ArcData object.
         */
        ArcData(node_t a, const int j, const Span<coord3d> X, const NodeNeighbours<K> &G)
        {
            __builtin_assume(j < 3);
            this->j = j;
            coord3d ap, am, ab, ac, ad, mp, db, bc;
            coord3d X_a = X[a], X_b = X[d_get(G.cubic_neighbours, j)];
            // Compute the arcs ab, ac, ad, bp, bm, ap, am, mp, bc and cd
            ab = (X_b - X_a);
            rabn = bond_length(ab);
            abh = rabn * ab;
            ac = (X[d_get(G.cubic_neighbours, (j + 1) % 3)] - X_a);
            racn = bond_length(ac);
            ach = racn * ac;
            rab = non_resciprocal_bond_length(ab);
            ad = (X[d_get(G.cubic_neighbours, (j + 2) % 3)] - X_a);
            radn = bond_length(ad);
            adh = radn * ad;
            db = (X_b - X[d_get(G.cubic_neighbours, (j + 2) % 3)]);
            rdbn = bond_length(db);
            dbh = rdbn * db;
            coord3d bp = (X[d_get(G.next_on_face, j)] - X_b);
            coord3d bm = (X[d_get(G.prev_on_face, j)] - X_b);

            rbpn = bond_length(bp);
            bph = bp * rbpn;
            rbmn = bond_length(bm);
            bmh = bm * rbmn;

            ap = bp + ab;
            rapn = bond_length(ap);
            aph = rapn * ap;
            am = bm + ab;
            ramn = bond_length(am);
            amh = ramn * am;
            mp = bp - bm;
            rmpn = bond_length(mp);
            mph = rmpn * mp;
            bc = ac - ab;
            rbcn = bond_length(bc);
            bch = rbcn * bc;
            cdh = unit_vector(ad - ac);

            // Compute inverses of some arcs, these are subject to be omitted if the equations are adapted appropriately with inversion of signs.
            bah = -abh;
            mbh = -bmh;
            pah = -aph;
            pbh = -bph;
        }

        // 3 FLOPs
        /**
         * @brief Compute the harmonic energy contribution from one parameter.
         * @param p0 The equillibrium value of the parameter.
         * @param p The current value of the parameter.
         */
        real_t harmonic_energy(const real_t p0, const real_t p) const
        {
            return (real_t)0.5 * (p - p0) * (p - p0);
        }
        // 4 FLOPs
        /**
         * @brief Compute the harmonic energy gradient contribution from one parameter.
         * @param p0 The equillibrium value of the parameter.
         * @param p The current value of the parameter.
         * @param gradp The gradient of the parameter with respect to the node position.
         */
        coord3d harmonic_energy_gradient(const real_t p0, const real_t p, const coord3d gradp) const
        {
            return (p - p0) * gradp;
        }

        mat3 harmonic_energy_hessian(const real_t p0, const real_t p, const coord3d grad_a, const coord3d grad_b, const mat3 &hessp) const
        {
            return hessp * (p - p0) + tensor_product(grad_a, grad_b);
        }

        // 1 FLOP
        /**
         * @brief Compute the bond length of the main arc ab or ac or ad. for j = 0, 1 or 2 respectively.
         * @return The bond length.
         */
        real_t bond() const { return rab; }
        // 5 FLOPs
        /**
         * @brief Compute the cosine of the angle between the main arc and the next arc, (ab,ac), (ac,ad), (ad,ab). For j = 0, 1 or 2 respectively.
         * @return The cosine of the angle.
         */
        real_t angle() const { return dot(abh, ach); }

        real_t normalized_angle_err() const { return sycl::acos(dot(abh, ach)); }

        // Returns outer angle m, used only diagnostically.
        real_t outer_angle_m() const { return -dot(abh, bmh); } // Compute outer angle. ab,bm

        // Returns outer angle p, used only diagnostically.
        real_t outer_angle_p() const { return -dot(abh, bph); } // Compute outer angle. ab,bp

        // Returns the inner dihedral angle for the current arc. Used here only for energy calculation,
        // otherwise embedded in dihedral computation because the planes and angles that make up the dihedral angle computation are required for derivative computation.
        // 50 FLOPs
        /**
         * @brief Compute the dihedral angle between the planes (abc,bcd), (acd,bcd) and (abd,bcd). For j = 0, 1 or 2 respectively.
         * @return The dihedral angle.
         */
        real_t dihedral() const
        {
            coord3d nabc, nbcd;
            real_t cos_b, cos_c, r_sin_b, r_sin_c;
            cos_b = dot(bah, bch);
            r_sin_b = sycl::rsqrt((real_t)1.0 - cos_b * cos_b);
            nabc = cross(bah, bch) * r_sin_b;
            cos_c = dot(-bch, cdh);
            r_sin_c = sycl::rsqrt((real_t)1.0 - cos_c * cos_c);
            nbcd = cross(-bch, cdh) * r_sin_c;
            return dot(nabc, nbcd);
        }
        // Returns the Outer-dihedral-a wrt. current arc, only accessed diagnostically (internal coordinate).
        /**
         * @brief Compute the dihedral angle between the planes $(b-a-b_m, a-b_m-b_p)$, $(c-a-c_m, a-c_m-c_p)$ and $(d-a-d_m, a-d_m-d_p)$. For j = 0, 1 or 2 respectively.
         * @return The dihedral angle.
         */
        real_t outer_dihedral_a() const
        {
            coord3d nbam_hat, namp_hat;
            real_t cos_a, cos_m, r_sin_a, r_sin_m;
            cos_a = dot(abh, amh);
            r_sin_a = sycl::rsqrt((real_t)1.0 - cos_a * cos_a);
            nbam_hat = cross(abh, amh) * r_sin_a;
            cos_m = dot(-amh, mph);
            r_sin_m = sycl::rsqrt((real_t)1.0 - cos_m * cos_m);
            namp_hat = cross(-amh, mph) * r_sin_m;
            real_t cos_beta = dot(nbam_hat, namp_hat); // Outer Dihedral angle bam, amp
            return cos_beta;
        }
        // Returns the Outer-dihedral-m wrt. current arc, only accessed diagnostically (internal coordinate).
        /**
         * @brief Compute the dihedral angle between the planes $(b-b_m-b_p, b_m-b_p-a)$, $(c-c_m-c_p, c_m-c_p-a)$ and $(d-d_m-d_p, d_m-d_p-a)$. For j = 0, 1 or 2 respectively.
         * @return The dihedral angle.
         */
        real_t outer_dihedral_m() const
        {
            coord3d nbmp_hat, nmpa_hat;
            real_t cos_m, cos_p, r_sin_m, r_sin_p;
            cos_m = dot(mbh, mph);
            r_sin_m = sycl::rsqrt((real_t)1.0 - cos_m * cos_m);
            nbmp_hat = cross(mbh, mph) * r_sin_m;
            cos_p = dot(-mph, pah);
            r_sin_p = sycl::rsqrt((real_t)1.0 - cos_p * cos_p);
            nmpa_hat = cross(-mph, pah) * r_sin_p;
            // Cosine to the outer dihedral angle constituted by the planes bmp and mpa
            real_t cos_beta = dot(nbmp_hat, nmpa_hat); // Outer dihedral angle bmp,mpa.
            return cos_beta;
        }
        // Returns the Outer-dihedral-p wrt. current arc, only accessed diagnostically (internal coordinate).
        /**
         * @brief Compute the dihedral angle between the planes $(b-b_p-a, b_p-a-b_m)$, $(c-c_p-a, c_p-a-c_m)$ and $(d-d_p-a, d_p-a-d_m)$. For j = 0, 1 or 2 respectively.
         * @return The dihedral angle.
         */
        real_t outer_dihedral_p() const
        {
            coord3d nbpa_hat, npam_hat;
            real_t cos_p, cos_a, r_sin_p, r_sin_a;
            cos_a = dot(aph, amh);
            r_sin_a = sycl::rsqrt((real_t)1.0 - cos_a * cos_a);
            npam_hat = cross(aph, amh) * r_sin_a;
            cos_p = dot(pbh, -aph);
            r_sin_p = sycl::rsqrt((real_t)1.0 - cos_p * cos_p);
            nbpa_hat = cross(pbh, -aph) * r_sin_p;
            real_t cos_beta = dot(nbpa_hat, npam_hat); // Outer dihedral angle bpa, pam.
            // Eq. 33 multiplied by harmonic term.
            return cos_beta;
        }

        // Chain rule terms for angle calculation
        // Computes gradient related to bending term. ~24 FLOPs
        /**
         * @brief Compute the gradient of the bending term.
         * @param c The constants for the threadIdx^th node.
         * @return The gradient of the bending term.
         */
        coord3d inner_angle_gradient(const Constants &c) const
        {
            real_t cos_angle = angle();                                                                       // Inner angle of arcs ab,ac.
            coord3d grad = cos_angle * (abh * rabn + ach * racn) - abh * racn - ach * rabn;                   // Derivative of inner angle: Eq. 21.
            return d_get(c.f_inner_angle, j) * harmonic_energy_gradient(d_get(c.angle0, j), cos_angle, grad); // Harmonic Energy Gradient: Eq. 21. multiplied by harmonic term.
        }

        mat3 bond_hessian_a(const Constants &c) const
        {
            coord3d grad_a = -abh;
            mat3 hessp = (identity3() - tensor_product(abh, abh)) * rabn;
            return d_get(c.f_bond, j) * harmonic_energy_hessian(d_get(c.r0, j), rab, grad_a, grad_a, hessp);
        }

        mat3 bond_hessian_b(const Constants &c) const
        {
            coord3d grad_a = -abh;
            coord3d grad_b = abh;
            mat3 hessp = (tensor_product(abh, abh) - identity3()) * rabn;
            return d_get(c.f_bond, j) * harmonic_energy_hessian(d_get(c.r0, j), rab, grad_a, grad_b, hessp);
        }

        mat3 inner_angle_hessian_a(const Constants &c) const
        {
            real_t cos_angle = angle(); // Inner angle of arcs ab,ac.
            coord3d grad_a = (abh * angle() - ach) * rabn + (ach * angle() - abh) * racn;

            mat3 GradG = tensor_product(abh, rabn * (abh * angle() - ach) + racn * (ach * angle() - abh)) + angle() * rabn * (tensor_product(abh, abh) - identity3()) - racn * (tensor_product(ach, ach) - identity3());

            mat3 GradF = tensor_product(ach, rabn * (abh * angle() - ach) + racn * (ach * angle() - abh)) + angle() * racn * (tensor_product(ach, ach) - identity3()) - rabn * (tensor_product(abh, abh) - identity3());

            mat3 P1 = tensor_product(abh * angle() - ach, abh * rabn * rabn);
            mat3 P2 = rabn * GradG;
            mat3 P3 = tensor_product(ach * angle() - abh, ach * racn * racn);
            mat3 P4 = racn * GradF;
            return c.f_inner_angle[j] * harmonic_energy_hessian(c.angle0[j], cos_angle, grad_a, grad_a, P1 + P2 + P3 + P4); // Harmonic Energy Hessian
        }

        mat3 inner_angle_hessian_b(const Constants &c) const
        {
            coord3d grad_a = (abh * angle() - ach) * rabn + (ach * angle() - abh) * racn;
            coord3d grad_b = (ach - abh * angle()) * rabn;
            mat3 G = tensor_product(abh, rabn * (ach - abh * angle())) + angle() * rabn * (identity3() - tensor_product(abh, abh));
            mat3 F = tensor_product(ach, rabn * (ach - abh * angle())) - rabn * (identity3() - tensor_product(abh, abh));
            mat3 P1 = tensor_product(abh * angle() - ach, -abh * rabn * rabn);

            mat3 P2 = rabn * G;
            mat3 P4 = racn * F;
            return d_get(c.f_inner_angle, j) * harmonic_energy_hessian(d_get(c.angle0, j), angle(), grad_a, grad_b, P1 + P2 + P4); // Harmonic Energy Hessian
        }

        mat3 inner_angle_hessian_c(const Constants &c) const
        {
            coord3d grad_a = (abh * angle() - ach) * rabn + (ach * angle() - abh) * racn;
            coord3d grad_c = (abh - angle() * ach) * racn;
            mat3 G = tensor_product(abh, racn * (abh - ach * angle())) - racn * (identity3() - tensor_product(ach, ach));
            mat3 F = tensor_product(ach, racn * (abh - ach * angle())) + angle() * racn * (identity3() - tensor_product(ach, ach));
            mat3 P2 = rabn * G;
            mat3 P3 = tensor_product(ach * angle() - abh, -ach * racn * racn);
            mat3 P4 = racn * F;
            return d_get(c.f_inner_angle, j) * harmonic_energy_hessian(d_get(c.angle0, j), angle(), grad_a, grad_c, P2 + P3 + P4); // Harmonic Energy Hessian
        }

        mat3 outer_angle_hessian_m_a(const Constants &c) const
        {
            real_t cost = dot(bah, bmh); // Compute outer angle. ab,bm
            mat3 gradba = rabn * (identity3() - tensor_product(bah, bah));
            coord3d grad_a = (bmh - bah * cost) * rabn;
            mat3 P1 = tensor_product(bmh - bah * cost, -bah * rabn * rabn);
            mat3 P2 = -rabn * (tensor_product(bah, grad_a) + cost * gradba);
            return d_get(c.f_outer_angle_m, j) * harmonic_energy_hessian(d_get(c.outer_angle_m0, j), cost, grad_a, grad_a, P1 + P2); // Harmonic Energy Hessian
        }

        mat3 outer_angle_hessian_m_b(const Constants &c) const
        {
            real_t cost = dot(bah, bmh); // Compute outer angle. ba,bm
            mat3 gradba = rabn * (tensor_product(bah, bah) - identity3());
            mat3 gradbm = rbmn * (tensor_product(bmh, bmh) - identity3());
            coord3d grad_b = rabn * (bah * cost - bmh) + rbmn * (bmh * cost - bah);
            coord3d grad_a = (bmh - bah * cost) * rabn;
            mat3 P1 = tensor_product(bmh - bah * cost, bah * rabn * rabn);
            mat3 P3 = rabn * (gradbm - (tensor_product(bah, grad_b) + cost * gradba));
            return d_get(c.f_outer_angle_m, j) * harmonic_energy_hessian(d_get(c.outer_angle_m0, j), cost, grad_a, grad_b, P1 + P3); // Harmonic Energy Hessian
        }

        mat3 outer_angle_hessian_m_m(const Constants &c) const
        {
            real_t cost = dot(bah, bmh); // Compute outer angle. ba,bm
            mat3 gradbm = rbmn * (identity3() - tensor_product(bmh, bmh));
            coord3d grad_a = (bmh - bah * cost) * rabn;
            coord3d grad_m = rbmn * (bah - bmh * cost);
            mat3 P1 = rabn * (gradbm - tensor_product(bah, grad_m));
            return d_get(c.f_outer_angle_m, j) * harmonic_energy_hessian(d_get(c.outer_angle_m0, j), cost, grad_a, grad_m, P1); // Harmonic Energy Hessian
        }

        mat3 outer_angle_hessian_p_a(const Constants &c) const
        {
            real_t cost = dot(bah, bph); // Compute outer angle. ba,bp
            mat3 gradba = rabn * (identity3() - tensor_product(bah, bah));
            coord3d grad_a = rabn * (bph - bah * cost);
            mat3 P1 = tensor_product(bph - bah * cost, -bah * rabn * rabn);
            mat3 P2 = -rabn * (tensor_product(bah, grad_a) + cost * gradba);
            return d_get(c.f_outer_angle_p, j) * harmonic_energy_hessian(d_get(c.outer_angle_p0, j), cost, grad_a, grad_a, P1 + P2); // Harmonic Energy Hessian
        }

        mat3 outer_angle_hessian_p_b(const Constants &c) const
        {
            real_t cost = dot(bah, bph); // Compute outer angle. ba,bp
            mat3 gradba = rabn * (tensor_product(bah, bah) - identity3());
            mat3 gradbp = rbpn * (tensor_product(bph, bph) - identity3());
            coord3d grad_b = rabn * (bah * cost - bph) + rbpn * (bph * cost - bah);
            coord3d grad_a = rabn * (bph - bah * cost);
            mat3 P1 = tensor_product(bph - bah * cost, bah * rabn * rabn);
            mat3 P3 = rabn * (gradbp - (tensor_product(bah, grad_b) + cost * gradba));
            return d_get(c.f_outer_angle_p, j) * harmonic_energy_hessian(d_get(c.outer_angle_p0, j), cost, grad_a, grad_b, P1 + P3); // Harmonic Energy Hessian
        }

        mat3 outer_angle_hessian_p_p(const Constants &c) const
        {
            real_t cost = dot(bah, bph); // Compute outer angle. ba,bp
            mat3 gradbp = rbpn * (identity3() - tensor_product(bph, bph));
            coord3d grad_a = rabn * (bph - bah * cost);
            coord3d grad_p = rbpn * (bah - bph * cost);
            mat3 P1 = rabn * (gradbp - tensor_product(bah, grad_p));
            return d_get(c.f_outer_angle_p, j) * harmonic_energy_hessian(d_get(c.outer_angle_p0, j), cost, grad_a, grad_p, P1); // Harmonic Energy Hessian
        }

        auto dihedral_hessian_terms(const Constants &c) const
        {
            coord3d cbh = -bch;
            real_t cost1 = dot(abh, cbh);
            real_t cost2 = dot(cbh, dbh);
            real_t sint1 = SQRT(1 - cost1 * cost1);
            real_t sint2 = SQRT(1 - cost2 * cost2);
            real_t cot1 = cost1 / sint1;
            real_t csc1 = real_t(1.) / sint1;
            real_t csc2 = real_t(1.) / sint2;
            coord3d nabc = cross(abh, cbh) * csc1;
            coord3d nbcd = cross(dbh, cbh) * csc2;
            real_t cosb = dot(nabc, nbcd);
            real_t Coeff = cosb * csc1 * rabn;
            coord3d F1 = abh * sint1;
            coord3d F2 = cross(cbh, nbcd) / cosb;
            coord3d F3 = cot1 * (abh * cost1 - cbh);
            coord3d F = F1 - F2 + F3;
            coord3d GradACosb = cosb * rabn * csc1 * (abh * sint1 - cross(cbh, nbcd) / cosb + cot1 * (abh * cost1 - cbh));
            return std::tuple{cbh, cost1, cost2, sint1, sint2, cot1, csc1, csc2, nabc, nbcd, cosb, Coeff, F, GradACosb};
        }

        // $\nabla_a(\nabla_a(\cos(\theta)))$
        mat3 dihedral_hessian_a(const Constants &c) const
        {
            auto [cbh, cost1, cost2, sint1, sint2, cot1, csc1, csc2, nabc, nbcd, cosb, Coeff, F, GradACosb] = dihedral_hessian_terms(c);
            coord3d GradARab = abh * rabn * rabn;
            mat3 GradAabh = (tensor_product(abh, abh) - identity3()) * rabn;
            coord3d GradASint1 = -(abh * cost1 - cbh) * cost1 * rabn * csc1;
            coord3d GradAcsc1 = -GradASint1 * csc1 * csc1;
            coord3d GradAcost1 = (abh * cost1 - cbh) * rabn;
            coord3d GradAcot1 = (sint1 * GradAcost1 - cost1 * GradASint1) * csc1 * csc1;
            coord3d GradACoeff = GradACosb * rabn * csc1 + cosb * (GradARab * csc1 + GradAcsc1 * rabn);
            mat3 GradAF1 = GradAabh * sint1 + tensor_product(abh, GradASint1);
            mat3 GradAF2 = tensor_product(cross(cbh, nbcd), -GradACosb / (cosb * cosb));
            mat3 GradAF3 = tensor_product(abh * cost1 - cbh, GradAcot1) + cot1 * (tensor_product(abh, GradAcost1) + cost1 * GradAabh);
            mat3 GradAF = GradAF1 - GradAF2 + GradAF3;
            mat3 GradAGradCosb = tensor_product(F, GradACoeff) + Coeff * GradAF;
            return d_get(c.f_inner_dihedral, j) * harmonic_energy_hessian(d_get(c.inner_dih0, j), cosb, GradACosb, GradACosb, GradAGradCosb); // Harmonic Energy Hessian
        }

        // $\nabla_b(\nabla_a(\cos(\beta)))$
        mat3 dihedral_hessian_b(const Constants &c) const
        {
            auto [cbh, cost1, cost2, sint1, sint2, cot1, csc1, csc2, nabc, nbcd, cosb, Coeff, F, GradACosb] = dihedral_hessian_terms(c);
            coord3d grad_b_sint1 = -((cbh - abh * cost1) * rabn + (abh - cbh * cost1) * rbcn) * cost1 * csc1;
            coord3d grad_b_sint2 = -((cbh - dbh * cost2) * rdbn + (dbh - cbh * cost2) * rbcn) * cost2 * csc2;
            coord3d grad_b_ab_cross_cb_dot_nbcd = (rbcn * (cross(nbcd, abh) - dot(nbcd, cross(abh, cbh)) * cbh) - rabn * (cross(nbcd, cbh) - dot(nbcd, cross(cbh, abh)) * abh));
            coord3d grad_b_db_cross_cb_dot_nabc = (rbcn * (cross(nabc, dbh) - dot(nabc, cross(dbh, cbh)) * cbh) - rdbn * (cross(nabc, cbh) - dot(nabc, cross(cbh, dbh)) * dbh));
            coord3d P1 = (grad_b_ab_cross_cb_dot_nbcd * sint1 - (dot(nbcd, cross(abh, cbh))) * grad_b_sint1) * csc1 * csc1;
            coord3d P2 = (grad_b_db_cross_cb_dot_nabc * sint2 - (dot(nabc, cross(dbh, cbh))) * grad_b_sint2) * csc2 * csc2;
            coord3d grad_b = P1 + P2;
            coord3d GradBRab = -abh * rabn * rabn;
            mat3 GradBabh = (identity3() - tensor_product(abh, abh)) * rabn;
            mat3 GradBcbh = (identity3() - tensor_product(cbh, cbh)) * rbcn;
            mat3 GradBdbh = (identity3() - tensor_product(dbh, dbh)) * rdbn;
            mat3 GradBnbcd = ((cross(dbh, GradBcbh) - cross(cbh, GradBdbh)) * sint2 - tensor_product(cross(dbh, cbh), grad_b_sint2)) * csc2 * csc2;
            coord3d GradBcsc1 = -grad_b_sint1 * csc1 * csc1;
            coord3d GradBcost1 = (cbh - abh * cost1) * rabn + (abh - cbh * cost1) * rbcn;
            coord3d GradBcot1 = (sint1 * GradBcost1 - cost1 * grad_b_sint1) * csc1 * csc1;
            coord3d GradBCoeff = grad_b * rabn * csc1 + cosb * (GradBRab * csc1 + GradBcsc1 * rabn);
            mat3 GradBF1 = GradBabh * sint1 + tensor_product(abh, grad_b_sint1);
            mat3 GradBF2 = tensor_product(cross(cbh, nbcd), -grad_b / (cosb * cosb)) + (cross(cbh, GradBnbcd) - cross(nbcd, GradBcbh)) / cosb;
            mat3 GradBF3 = tensor_product(abh * cost1 - cbh, GradBcot1) + cot1 * (tensor_product(abh, GradBcost1) + cost1 * GradBabh - GradBcbh);
            mat3 GradBF = GradBF1 - GradBF2 + GradBF3;
            mat3 GradBGradCosb = tensor_product(F, GradBCoeff) + Coeff * GradBF;
            return d_get(c.f_inner_dihedral, j) * harmonic_energy_hessian(d_get(c.inner_dih0, j), cosb, GradACosb, grad_b, GradBGradCosb); // Harmonic Energy Hessian
        }

        // $\nabla_c(\nabla_a(\cos(\theta)))$
        mat3 dihedral_hessian_c(const Constants &c) const
        {
            auto [cbh, cost1, cost2, sint1, sint2, cot1, csc1, csc2, nabc, nbcd, cosb, Coeff, F, GradACosb] = dihedral_hessian_terms(c);
            coord3d grad_c_sint1 = -(cbh * cost1 - abh) * cost1 * csc1 * rbcn;
            coord3d grad_c_sint2 = -(cbh * cost2 - dbh) * cost2 * csc2 * rbcn;
            coord3d grad_c_ab_cross_cb_dot_nabc = rbcn * (dot(nabc, cross(dbh, cbh)) * cbh - cross(nabc, dbh));
            coord3d grad_c_db_cross_cb_dot_nbcd = rbcn * (dot(nbcd, cross(abh, cbh)) * cbh - cross(nbcd, abh));
            coord3d P1 = (grad_c_ab_cross_cb_dot_nabc * sint2 - (dot(nabc, cross(dbh, cbh))) * grad_c_sint2) * csc2 * csc2;
            coord3d P2 = (grad_c_db_cross_cb_dot_nbcd * sint1 - (dot(nbcd, cross(abh, cbh))) * grad_c_sint1) * csc1 * csc1;
            coord3d grad_c = P1 + P2;

            mat3 GradCcbh = (tensor_product(cbh, cbh) - identity3()) * rbcn;
            coord3d GradCcsc1 = -grad_c_sint1 * csc1 * csc1;
            coord3d GradCcost1 = (cbh * cost1 - abh) * rbcn;
            coord3d GradCcot1 = (sint1 * GradCcost1 - cost1 * grad_c_sint1) * csc1 * csc1;
            mat3 GradCnbcd = (cross(dbh, GradCcbh) * sint2 - tensor_product(cross(dbh, cbh), grad_c_sint2)) * csc2 * csc2;
            coord3d GradCCoeff = grad_c * rabn * csc1 + cosb * (GradCcsc1 * rabn);
            mat3 GradCF1 = tensor_product(abh, grad_c_sint1);
            mat3 GradCF2 = tensor_product(cross(cbh, nbcd), -grad_c / (cosb * cosb)) + (cross(cbh, GradCnbcd) - cross(nbcd, GradCcbh)) / cosb;
            mat3 GradCF3 = tensor_product(abh * cost1 - cbh, GradCcot1) + cot1 * (tensor_product(abh, GradCcost1) - GradCcbh);
            mat3 GradCF = GradCF1 - GradCF2 + GradCF3;
            mat3 GradCGradCosb = tensor_product(F, GradCCoeff) + Coeff * GradCF;
            return d_get(c.f_inner_dihedral, j) * harmonic_energy_hessian(d_get(c.inner_dih0, j), cosb, GradACosb, grad_c, GradCGradCosb); // Harmonic Energy Hessian
        }

        // $\nabla_d(\nabla_a(\cos(\theta)))$
        mat3 dihedral_hessian_d(const Constants &c) const
        {
            auto [cbh, cost1, cost2, sint1, sint2, cot1, csc1, csc2, nabc, nbcd, cosb, Coeff, F, GradACosb] = dihedral_hessian_terms(c);
            coord3d GradDSint2 = -(dbh * cost2 - cbh) * cost2 * csc2 * rdbn;
            coord3d GradDDbCrossCbDotNabc = -rdbn * (dot(nabc, cross(cbh, dbh)) * dbh - cross(nabc, cbh));
            coord3d grad_d = (GradDDbCrossCbDotNabc * sint2 - (dot(nabc, cross(dbh, cbh))) * GradDSint2) * csc2 * csc2;

            mat3 GradDdbh = (tensor_product(dbh, dbh) - identity3()) * rdbn;
            mat3 GradDnbcd = (-cross(cbh, GradDdbh) * sint2 - tensor_product(cross(dbh, cbh), GradDSint2)) * csc2 * csc2;
            coord3d GradDCoeff = grad_d * rabn * csc1;
            mat3 GradDF2 = tensor_product(cross(cbh, nbcd), -grad_d / (cosb * cosb)) + cross(cbh, GradDnbcd) / cosb;
            mat3 GradDF = -GradDF2;
            mat3 GradDGradCosb = tensor_product(F, GradDCoeff) + Coeff * GradDF;
            return d_get(c.f_inner_dihedral, j) * harmonic_energy_hessian(d_get(c.inner_dih0, j), cosb, GradACosb, grad_d, GradDGradCosb); // Harmonic Energy Hessian
        }

        auto outer_dihedral_hessian_a_terms(const Constants &c) const
        {
            coord3d mah = -amh;
            real_t cosa = dot(abh, amh);
            real_t cosm = dot(mah, mph);
            real_t sina = sqrt(1 - cosa * cosa);
            real_t sinm = sqrt(1 - cosm * cosm);
            real_t cota = cosa / sina;
            real_t cotm = cosm / sinm;
            real_t csca = real_t(1.) / sina;
            real_t cscm = real_t(1.) / sinm;
            coord3d nbam = cross(abh, amh) * csca;
            coord3d namp = cross(mah, mph) * cscm;
            real_t cosb = dot(nbam, namp);
            coord3d F1 = abh * cosb;
            coord3d F2 = cross(amh, namp) * csca;
            coord3d F3 = amh * cosb;
            coord3d F4 = cross(namp, abh) * csca;
            coord3d G1 = abh * cosa * rabn;
            coord3d G2 = amh * rabn;
            coord3d G3 = amh * cosa * ramn;
            coord3d G4 = abh * ramn;
            coord3d H1 = cross(mph, nbam);
            coord3d H2 = mah * cosb * sinm;
            coord3d H3 = cotm * cosb * (mph - mah * cosm);
            real_t C1 = cota * cosb * csca;
            real_t C2 = ramn * cscm;
            coord3d GradAcosb = (F1 - F2) * rabn + (F3 - F4) * ramn + C1 * (G1 - G2 + G3 - G4) + C2 * (H1 - H2 + H3);
            return std::tuple(cosa, cosm, sina, sinm, cota, cotm, csca, cscm, nbam, namp, cosb, F1, F2, F3, F4, G1, G2, G3, G4, H1, H2, H3, C1, C2, GradAcosb);
        }

        mat3 outer_dihedral_hessian_a_a(const Constants &c) const
        {
            auto [cosa, cosm, sina, sinm, cota, cotm, csca, cscm, nbam, namp, cosb, F1, F2, F3, F4, G1, G2, G3, G4, H1, H2, H3, C1, C2, GradAcosb] = outer_dihedral_hessian_a_terms(c);
            coord3d mah = -amh;

            coord3d GradAcosa = (amh * cosa - abh) * ramn + (abh * cosa - amh) * rabn;
            coord3d GradAsina = -cosa * csca * GradAcosa;
            coord3d GradAcota = (sina * GradAcosa - cosa * GradAsina) * csca * csca;
            coord3d GradAcsca = -GradAsina * csca * csca;
            coord3d GradAcosm = (mph - mah * cosm) * ramn;
            coord3d GradAsinm = -cosm * cscm * GradAcosm;
            coord3d GradAcotm = (sinm * GradAcosm - cosm * GradAsinm) * cscm * cscm;
            coord3d GradAcscm = -GradAsinm * cscm * cscm;

            mat3 GradAab = (tensor_product(abh, abh) - identity3()) * rabn;
            coord3d GradArabn = abh * rabn * rabn;
            mat3 GradAam = (tensor_product(amh, amh) - identity3()) * ramn;
            coord3d GradAramn = amh * ramn * ramn;
            mat3 GradAma = (identity3() - tensor_product(mah, mah)) * ramn;
            mat3 GradAnbam = ((cross(abh, GradAam) - cross(amh, GradAab)) * sina - tensor_product(cross(abh, amh), GradAsina)) * csca * csca;
            mat3 GradAnamp = ((-cross(mph, GradAma)) * sinm - tensor_product(cross(mah, mph), GradAsinm)) * cscm * cscm;
            mat3 GradAF1 = tensor_product(abh, GradAcosb) + cosb * GradAab;
            mat3 GradAF2 = ((cross(amh, GradAnamp) - cross(namp, GradAam)) * sina - tensor_product(cross(amh, namp), GradAsina)) * csca * csca;
            mat3 GradAF3 = tensor_product(amh, GradAcosb) + cosb * GradAam;
            mat3 GradAF4 = ((cross(namp, GradAab) - cross(abh, GradAnamp)) * sina - tensor_product(cross(namp, abh), GradAsina)) * csca * csca;

            mat3 GradAG1 = tensor_product(abh, (GradAcosa * rabn + GradArabn * cosa)) + cosa * rabn * GradAab;
            mat3 GradAG2 = tensor_product(amh, GradArabn) + GradAam * rabn;
            mat3 GradAG3 = tensor_product(amh, (GradAcosa * ramn + GradAramn * cosa)) + cosa * ramn * GradAam;
            mat3 GradAG4 = tensor_product(abh, GradAramn) + GradAab * ramn;

            mat3 GradAH1 = cross(mph, GradAnbam);
            mat3 GradAH2 = tensor_product(mah, (cosb * GradAsinm + GradAcosb * sinm)) + cosb * sinm * GradAma;
            mat3 GradAH3 = tensor_product(mph - mah * cosm, (GradAcotm * cosb + GradAcosb * cotm)) + cotm * cosb * (-(GradAma * cosm + tensor_product(mah, GradAcosm)));

            coord3d GradAC1 = GradAcota * cosb * csca + cota * (GradAcosb * csca + cosb * GradAcsca);
            coord3d GradAC2 = GradAramn * cscm + GradAcscm * ramn;

            mat3 GradAF = rabn * (GradAF1 - GradAF2) + tensor_product(F1 - F2, GradArabn) + ramn * (GradAF3 - GradAF4) + tensor_product(F3 - F4, GradAramn);
            mat3 GradAG = C1 * (GradAG1 - GradAG2 + GradAG3 - GradAG4) + tensor_product(G1 - G2 + G3 - G4, GradAC1);
            mat3 GradAH = C2 * (GradAH1 - GradAH2 + GradAH3) + tensor_product(H1 - H2 + H3, GradAC2);

            mat3 GradGradAcosb = GradAF + GradAG + GradAH;
            return d_get(c.f_outer_dihedral, j) * harmonic_energy_hessian(d_get(c.outer_dih0_a, j), cosb, GradAcosb, GradAcosb, GradGradAcosb); // Harmonic Energy Hessian
        }

        mat3 outer_dihedral_hessian_a_b(const Constants &c) const
        {
            auto [cosa, cosm, sina, sinm, cota, cotm, csca, cscm, nbam, namp, cosb, F1, F2, F3, F4, G1, G2, G3, G4, H1, H2, H3, C1, C2, GradAcosb] = outer_dihedral_hessian_a_terms(c);
            coord3d mah = -amh;

            mat3 GradBab = (identity3() - tensor_product(abh, abh)) * rabn;
            coord3d GradBrabn = -abh * rabn * rabn;
            coord3d GradBcosa = (amh - cosa * abh) * rabn;
            coord3d GradBsina = -cosa / sina * GradBcosa;
            coord3d GradBcota = (sina * GradBcosa - cosa * GradBsina) * csca * csca;
            coord3d GradBcsca = -GradBsina * csca * csca;
            mat3 GradBnbam = ((-cross(amh, GradBab)) * sina - tensor_product(cross(abh, amh), GradBsina)) * csca * csca;
            coord3d GradBcosb = (-(cross(namp, amh) - dot(namp, cross(amh, abh)) * abh) * sina * rabn - dot(namp, cross(abh, amh)) * GradBsina) * csca * csca;

            mat3 GradBF1 = tensor_product(abh, GradBcosb) + cosb * GradBab;
            mat3 GradBF2 = (-tensor_product(cross(amh, namp), GradBsina)) * csca * csca;
            mat3 GradBF3 = tensor_product(amh, GradBcosb);
            mat3 GradBF4 = ((cross(namp, GradBab)) * sina - tensor_product(cross(namp, abh), GradBsina)) * csca * csca;

            mat3 GradBG1 = tensor_product(abh, (GradBcosa * rabn + GradBrabn * cosa)) + cosa * rabn * GradBab;
            mat3 GradBG2 = tensor_product(amh, GradBrabn);
            mat3 GradBG3 = tensor_product(amh, (GradBcosa * ramn));
            mat3 GradBG4 = GradBab * ramn;

            mat3 GradBH1 = cross(mph, GradBnbam);
            mat3 GradBH2 = tensor_product(mah, GradBcosb * sinm);
            mat3 GradBH3 = tensor_product(mph - mah * cosm, (GradBcosb * cotm));

            coord3d GradBC1 = GradBcota * cosb * csca + cota * (GradBcosb * csca + cosb * GradBcsca);

            mat3 GradBF = rabn * (GradBF1 - GradBF2) + tensor_product(F1 - F2, GradBrabn) + ramn * (GradBF3 - GradBF4);
            mat3 GradBG = C1 * (GradBG1 - GradBG2 + GradBG3 - GradBG4) + tensor_product(G1 - G2 + G3 - G4, GradBC1);
            mat3 GradBH = C2 * (GradBH1 - GradBH2 + GradBH3);

            mat3 GradGradBcosb = GradBF + GradBG + GradBH;
            return d_get(c.f_outer_dihedral, j) * harmonic_energy_hessian(d_get(c.outer_dih0_a, j), cosb, GradAcosb, GradBcosb, GradGradBcosb); // Harmonic Energy Hessian
        }

        mat3 outer_dihedral_hessian_a_m(const Constants &c) const
        {
            auto [cosa, cosm, sina, sinm, cota, cotm, csca, cscm, nbam, namp, cosb, F1, F2, F3, F4, G1, G2, G3, G4, H1, H2, H3, C1, C2, GradAcosb] = outer_dihedral_hessian_a_terms(c);
            coord3d mah = -amh;

            coord3d GradMcosa = (abh - cosa * amh) * ramn;
            coord3d GradMcosm = (mah * cosm - mph) * ramn + (mph * cosm - mah) * rmpn;
            coord3d GradMsina = -cosa / sina * GradMcosa;
            coord3d GradMsinm = -cosm / sinm * GradMcosm;

            mat3 GradMam = (identity3() - tensor_product(amh, amh)) * ramn;
            coord3d GradMramn = -amh * ramn * ramn;
            mat3 GradMma = (tensor_product(amh, amh) - identity3()) * ramn;
            mat3 GradMmp = (tensor_product(mph, mph) - identity3()) * rmpn;
            coord3d GradMcota = (sina * GradMcosa - cosa * GradMsina) * csca * csca;
            coord3d GradMcotm = (sinm * GradMcosm - cosm * GradMsinm) * cscm * cscm;
            coord3d GradMcsca = -GradMsina * csca * csca;
            coord3d GradMcscm = -GradMsinm * cscm * cscm;
            mat3 GradMnbam = ((cross(abh, GradMam)) * sina - tensor_product(cross(abh, amh), GradMsina)) * csca * csca;
            mat3 GradMnamp = ((cross(mah, GradMmp) - cross(mph, GradMma)) * sinm - tensor_product(cross(mah, mph), GradMsinm)) * cscm * cscm;

            coord3d cosbP1 = (((cross(namp, abh) - dot(namp, cross(abh, amh)) * amh) * sina * ramn - dot(namp, cross(abh, amh)) * GradMsina) * csca * csca);
            coord3d cosbP2 = (((rmpn * cscm) * (dot(nbam, cross(mah, mph)) * mph - cross(nbam, mah)) - (cscm * ramn) * (dot(nbam, cross(mph, mah)) * mah - cross(nbam, mph))) - dot(nbam, cross(mah, mph)) * GradMsinm * cscm * cscm);
            coord3d GradMcosb = cosbP1 + cosbP2;
            mat3 GradMF1 = tensor_product(abh, GradMcosb);
            mat3 GradMF2 = ((cross(amh, GradMnamp) - cross(namp, GradMam)) * sina - tensor_product(cross(amh, namp), GradMsina)) * csca * csca;
            mat3 GradMF3 = tensor_product(amh, GradMcosb) + cosb * GradMam;
            mat3 GradMF4 = ((-cross(abh, GradMnamp)) * sina - tensor_product(cross(namp, abh), GradMsina)) * csca * csca;

            mat3 GradMG1 = tensor_product(abh, GradMcosa * rabn);
            mat3 GradMG2 = GradMam * rabn;
            mat3 GradMG3 = tensor_product(amh, GradMcosa * ramn + GradMramn * cosa) + cosa * ramn * GradMam;
            mat3 GradMG4 = tensor_product(abh, GradMramn);

            mat3 GradMH1 = cross(mph, GradMnbam) - cross(nbam, GradMmp);
            mat3 GradMH2 = tensor_product(mah, (cosb * GradMsinm + GradMcosb * sinm)) + cosb * sinm * GradMma;
            mat3 GradMH3 = tensor_product(mph - mah * cosm, (GradMcotm * cosb + GradMcosb * cotm)) + cotm * cosb * (GradMmp - (GradMma * cosm + tensor_product(mah, GradMcosm)));

            coord3d GradMC1 = GradMcota * cosb * csca + cota * (GradMcosb * csca + cosb * GradMcsca);
            coord3d GradMC2 = GradMramn * cscm + GradMcscm * ramn;

            mat3 GradMF = rabn * (GradMF1 - GradMF2) + ramn * (GradMF3 - GradMF4) + tensor_product(F3 - F4, GradMramn);
            mat3 GradMG = C1 * (GradMG1 - GradMG2 + GradMG3 - GradMG4) + tensor_product(G1 - G2 + G3 - G4, GradMC1);
            mat3 GradMH = C2 * (GradMH1 - GradMH2 + GradMH3) + tensor_product(H1 - H2 + H3, GradMC2);

            mat3 GradGradMcosb = GradMF + GradMG + GradMH;
            return d_get(c.f_outer_dihedral, j) * harmonic_energy_hessian(d_get(c.outer_dih0_a, j), cosb, GradAcosb, GradMcosb, GradGradMcosb); // Harmonic Energy Hessian
        }

        mat3 outer_dihedral_hessian_a_p(const Constants &c) const
        {
            auto [cosa, cosm, sina, sinm, cota, cotm, csca, cscm, nbam, namp, cosb, F1, F2, F3, F4, G1, G2, G3, G4, H1, H2, H3, C1, C2, GradAcosb] = outer_dihedral_hessian_a_terms(c);
            coord3d mah = -amh;

            coord3d GradPcosm = (mah - mph * cosm) * rmpn;
            coord3d GradPsinm = -cosm * cscm * GradPcosm;

            mat3 GradPmp = (identity3() - tensor_product(mph, mph)) * rmpn;
            coord3d GradPcotm = (sinm * GradPcosm - cosm * GradPsinm) * cscm * cscm;
            coord3d GradPcscm = -GradPsinm * cscm * cscm;
            mat3 GradPnbam = -tensor_product(cross(abh, amh), GradAcosb) * csca * csca;
            mat3 GradPnamp = (cross(mah, GradPmp) * sinm - tensor_product(cross(mah, mph), GradPsinm)) * cscm * cscm;
            coord3d GradPcosb = (sinm * rmpn * (cross(nbam, mah) - dot(nbam, cross(mah, mph)) * mph) - dot(nbam, cross(mah, mph)) * GradPsinm) * cscm * cscm;

            mat3 GradPF1 = tensor_product(abh, GradPcosb);
            mat3 GradPF2 = ((cross(amh, GradPnamp)) * sina) * csca * csca;
            mat3 GradPF3 = tensor_product(amh, GradPcosb);
            mat3 GradPF4 = ((-cross(abh, GradPnamp)) * sina) * csca * csca;

            mat3 GradPH1 = -cross(nbam, GradPmp);
            mat3 GradPH2 = tensor_product(mah, (cosb * GradPsinm + GradPcosb * sinm));
            mat3 GradPH3 = tensor_product(mph - mah * cosm, (GradPcotm * cosb + GradPcosb * cotm)) + cotm * cosb * (GradPmp - (tensor_product(mah, GradPcosm)));

            coord3d GradPC1 = cota * (GradPcosb * csca);
            coord3d GradPC2 = GradPcscm * ramn;

            mat3 GradPF = (GradPF1 - GradPF2) * rabn + (GradPF3 - GradPF4) * ramn;
            mat3 GradPG = tensor_product(G1 - G2 + G3 - G4, GradPC1);
            mat3 GradPH = C2 * (GradPH1 - GradPH2 + GradPH3) + tensor_product(H1 - H2 + H3, GradPC2);

            mat3 GradGradPcosb = GradPF + GradPG + GradPH;
            return d_get(c.f_outer_dihedral, j) * harmonic_energy_hessian(d_get(c.outer_dih0_a, j), cosb, GradAcosb, GradPcosb, GradGradPcosb); // Harmonic Energy Hessian
        }

        auto outer_dihedral_hessian_m_terms() const
        {
            coord3d pmh = -mph;
            real_t cosm = dot(mbh, mph);
            real_t cosp = dot(pmh, pah);
            real_t sinm = sqrt(1 - cosm * cosm);
            real_t sinp = sqrt(1 - cosp * cosp);
            real_t cscm = 1 / sinm;
            real_t cscp = 1 / sinp;
            real_t cotp = cosp / sinp;
            real_t cotm = cosm / sinm;
            coord3d nbmp = cross(mbh, mph) * cscm;
            coord3d nmpa = cross(pmh, pah) * cscp;
            real_t cosb = dot(nbmp, nmpa);
            coord3d F = cross(nbmp, pmh) * rapn * cscp;
            coord3d G = pah * cosb * rapn;
            real_t K1 = cotp * cosb;
            real_t K2 = rapn * cscp;
            real_t K_ = K1 * K2;
            coord3d H = K_ * (pmh - pah * cosp);
            coord3d GradAcosb = F - G + H;
            return std::tuple(cosm, cosp, sinm, sinp, cscm, cscp, cotp, cotm, nbmp, nmpa, cosb, F, G, H, K1, K2, K_, GradAcosb);
        }

        mat3 outer_dihedral_hessian_m_a(const Constants &c) const
        {
            auto [cosm, cosp, sinm, sinp, cscm, cscp, cotp, cotm, nbmp, nmpa, cosb, F, G, H, K1, K2, K_, GradAcosb] = outer_dihedral_hessian_m_terms();
            coord3d pmh = -mph;

            coord3d GradAcosp = (pmh - pah * cosp) * rapn;
            coord3d GradAsinp = -cosp * cscp * GradAcosp;
            coord3d GradAcscp = -GradAsinp * cscp * cscp;
            coord3d GradAcotp = (sinp * GradAcosp - cosp * GradAsinp) * cscp * cscp;

            mat3 GradApah = (identity3() - tensor_product(pah, pah)) * rapn;
            coord3d GradArpan = -pah * rapn * rapn;
            mat3 GradAF = tensor_product(cross(nbmp, pmh), GradArpan * cscp + GradAcscp * rapn);
            mat3 GradAG = tensor_product(pah, GradAcosb * rapn + GradArpan * cosb) + GradApah * cosb * rapn;
            coord3d GradAK1 = GradAcotp * cosb + cotp * GradAcosb;
            coord3d GradAK2 = GradArpan * cscp + GradAcscp * rapn;
            coord3d GradAK = GradAK1 * K2 + K1 * GradAK2;
            mat3 GradAH = tensor_product(pmh - pah * cosp, GradAK) + K_ * (-GradApah * cosp - tensor_product(pah, GradAcosp));

            mat3 GradGradAcosb = GradAF - GradAG + GradAH;
            return d_get(c.f_outer_dihedral, j) * harmonic_energy_hessian(d_get(c.outer_dih0_m, j), cosb, GradAcosb, GradAcosb, GradGradAcosb); // Harmonic Energy Hessian
        }

        mat3 outer_dihedral_hessian_m_b(const Constants &c) const
        {
            auto [cosm, cosp, sinm, sinp, cscm, cscp, cotp, cotm, nbmp, nmpa, cosb, F, G, H, K1, K2, K_, GradAcosb] = outer_dihedral_hessian_m_terms();
            coord3d pmh = -mph;

            coord3d GradBcosm = (mph - mbh * cosm) * rbmn;
            coord3d GradBsinm = -cosm / sinm * GradBcosm;
            coord3d GradBcscm = -GradBsinm * cscm * cscm;

            coord3d GradBcosb = -(cross(nmpa, mph) - dot(cross(nmpa, mph), mbh) * mbh) * rbmn * cscm - dot(nmpa, cross(mbh, mph)) * GradBsinm * cscm * cscm;

            mat3 GradBmbh = (identity3() - tensor_product(mbh, mbh)) * rbmn;
            mat3 GradBnbmp = (-cross(mph, GradBmbh) * sinm - tensor_product(cross(mbh, mph), GradBsinm)) * cscm * cscm;
            mat3 GradBF = -cross(pmh, GradBnbmp) * rapn * cscp;
            mat3 GradBG = tensor_product(pah, GradBcosb * rapn);
            coord3d GradBK1 = cotp * GradBcosb;
            coord3d GradBK = GradBK1 * K2;
            mat3 GradBH = tensor_product(pmh - pah * cosp, GradBK);
            mat3 GradGradBcosb = GradBF - GradBG + GradBH;
            return d_get(c.f_outer_dihedral, j) * harmonic_energy_hessian(d_get(c.outer_dih0_m, j), cosb, GradAcosb, GradBcosb, GradGradBcosb); // Harmonic Energy Hessian
        }

        mat3 outer_dihedral_hessian_m_m(const Constants &c) const
        {
            auto [cosm, cosp, sinm, sinp, cscm, cscp, cotp, cotm, nbmp, nmpa, cosb, F, G, H, K1, K2, K_, GradAcosb] = outer_dihedral_hessian_m_terms();
            coord3d pmh = -mph;

            coord3d GradMcosm = (mbh * cosm - mph) * rbmn + (mph * cosm - mbh) * rmpn;
            coord3d GradMsinm = -cosm / sinm * GradMcosm;
            coord3d GradMcscm = -GradMsinm * cscm * cscm;
            coord3d GradMcosp = (pah - pmh * cosp) * rmpn;
            coord3d GradMsinp = -cosp * cscp * GradMcosp;
            coord3d GradMcscp = -GradMsinp * cscp * cscp;
            coord3d GradMcotp = (sinp * GradMcosp - cosp * GradMsinp) * cscp * cscp;

            mat3 GradMmph = (tensor_product(mph, mph) - identity3()) * rmpn;
            mat3 GradMpmh = (identity3() - tensor_product(pmh, pmh)) * rmpn;
            mat3 GradMmbh = (tensor_product(mbh, mbh) - identity3()) * rbmn;

            mat3 GradMnbmp = ((cross(mbh, GradMmph) - cross(mph, GradMmbh)) * sinm - tensor_product(cross(mbh, mph), GradMsinm)) * cscm * cscm;
            coord3d GradMcosbP1 = (((dot(cross(nmpa, mbh), mph) * mph - cross(nmpa, mbh)) * rmpn - (dot(cross(nmpa, mph), mbh) * mbh - cross(nmpa, mph)) * rbmn) * sinm - dot(nmpa, cross(mbh, mph)) * GradMsinm) * cscm * cscm;
            coord3d GradMcosbP2 = ((-(cross(nbmp, pah) - dot(cross(nbmp, pah), pmh) * pmh) * rmpn) * sinp - dot(nbmp, cross(pmh, pah)) * GradMsinp) * cscp * cscp;
            coord3d GradMcosb = GradMcosbP1 + GradMcosbP2;

            mat3 GradMF = (cross(nbmp, GradMpmh) - cross(pmh, GradMnbmp)) * rapn * cscp + tensor_product(cross(nbmp, pmh), GradMcscp * rapn);
            mat3 GradMG = tensor_product(pah, GradMcosb * rapn);
            coord3d GradMK1 = GradMcotp * cosb + cotp * GradMcosb;
            coord3d GradMK2 = GradMcscp * rapn;
            coord3d GradMK = GradMK1 * K2 + K1 * GradMK2;
            mat3 GradMH = tensor_product((pmh - pah * cosp), GradMK) + K_ * (GradMpmh - tensor_product(pah, GradMcosp));
            mat3 GradGradMcosb = GradMF - GradMG + GradMH;
            return d_get(c.f_outer_dihedral, j) * harmonic_energy_hessian(d_get(c.outer_dih0_m, j), cosb, GradAcosb, GradMcosb, GradGradMcosb); // Harmonic Energy Hessian
        }

        mat3 outer_dihedral_hessian_m_p(const Constants &c) const
        {
            auto [cosm, cosp, sinm, sinp, cscm, cscp, cotp, cotm, nbmp, nmpa, cosb, F, G, H, K1, K2, K_, GradAcosb] = outer_dihedral_hessian_m_terms();
            coord3d pmh = -mph;

            coord3d GradPcosm = (mbh - mph * cosm) * rmpn;
            coord3d GradPsinm = -cosm / sinm * GradPcosm;
            coord3d GradPcscm = -GradPsinm * cscm * cscm;
            coord3d GradPcosp = (pmh * cosp - pah) * rmpn + (pah * cosp - pmh) * rapn;
            coord3d GradPsinp = -cosp * cscp * GradPcosp;
            coord3d GradPcscp = -GradPsinp * cscp * cscp;
            coord3d GradPcotp = (sinp * GradPcosp - cosp * GradPsinp) * cscp * cscp;

            mat3 GradPmph = (identity3() - tensor_product(mph, mph)) * rmpn;
            mat3 GradPpmh = (tensor_product(pmh, pmh) - identity3()) * rmpn;
            mat3 GradPpah = (tensor_product(pah, pah) - identity3()) * rapn;

            coord3d GradPrpan = pah * rapn * rapn;

            mat3 GradPnbmp = ((cross(mbh, GradPmph)) * sinm - tensor_product(cross(mbh, mph), GradPsinm)) * cscm * cscm;
            coord3d GradPcosbP1 = (((cross(nmpa, mbh) - dot(cross(nmpa, mbh), mph) * mph) * rmpn) * sinm - dot(nmpa, cross(mbh, mph)) * GradPsinm) * cscm * cscm;
            coord3d GradPcosbP2 = (((dot(cross(nbmp, pmh), pah) * pah - cross(nbmp, pmh)) * rapn - (dot(cross(nbmp, pah), pmh) * pmh - cross(nbmp, pah)) * rmpn) * sinp - dot(nbmp, cross(pmh, pah)) * GradPsinp) * cscp * cscp;
            coord3d GradPcosb = GradPcosbP1 + GradPcosbP2;

            mat3 GradPF = (cross(nbmp, GradPpmh) - cross(pmh, GradPnbmp)) * rapn * cscp + tensor_product(cross(nbmp, pmh), GradPcscp * rapn + GradPrpan * cscp);
            mat3 GradPG = tensor_product(pah, GradPcosb * rapn + GradPrpan * cosb) + GradPpah * rapn * cosb;
            coord3d GradPK1 = GradPcotp * cosb + cotp * GradPcosb;
            coord3d GradPK2 = GradPrpan * cscp + GradPcscp * rapn;
            coord3d GradPK = GradPK1 * K2 + K1 * GradPK2;
            mat3 GradPH = tensor_product((pmh - pah * cosp), GradPK) + K_ * (GradPpmh - tensor_product(pah, GradPcosp) - GradPpah * cosp);

            mat3 GradGradPcosb = GradPF - GradPG + GradPH;
            return d_get(c.f_outer_dihedral, j) * harmonic_energy_hessian(d_get(c.outer_dih0_m, j), cosb, GradAcosb, GradPcosb, GradGradPcosb); // Harmonic Energy Hessian
        }

        auto outer_dihedral_hessian_p_terms() const
        {
            coord3d pah = -aph;
            real_t cosa = dot(aph, amh);
            real_t cosp = dot(pbh, pah);
            real_t sina = sqrt(1 - cosa * cosa);
            real_t sinp = sqrt(1 - cosp * cosp);
            real_t csca = 1 / sina;
            real_t cscp = 1 / sinp;
            real_t cotp = cosp / sinp;
            real_t cota = cosa / sina;
            coord3d nbpa = cross(pbh, pah) * cscp;
            coord3d npam = cross(aph, amh) * csca;
            real_t cosb = dot(nbpa, npam);

            real_t rpan = rapn;
            real_t C1 = rpan * cscp;
            real_t C2 = cota * cosb * csca;
            coord3d F1 = cross(npam, pbh);
            coord3d F2 = pah * cosb * sinp;
            coord3d F3 = cotp * cosb * (pbh - pah * cosp);
            coord3d G1 = aph * cosb;
            coord3d G2 = cross(amh, nbpa) * csca;
            coord3d G3 = amh * cosb;
            coord3d G4 = cross(nbpa, aph) * csca;
            coord3d H1 = aph * cosa * rapn;
            coord3d H2 = amh * rapn;
            coord3d H3 = amh * cosa * ramn;
            coord3d H4 = aph * ramn;

            coord3d GradAcosb = C1 * (F1 - F2 + F3) + rapn * (G1 - G2) + ramn * (G3 - G4) + C2 * (H1 - H2 + H3 - H4);
            return std::tuple(cosa, cosp, sina, sinp, csca, cscp, cotp, cota, nbpa, npam, cosb, C1, C2, F1, F2, F3, G1, G2, G3, G4, H1, H2, H3, H4, GradAcosb);
        }

        mat3 outer_dihedral_hessian_p_a(const Constants &c) const
        {
            auto [cosa, cosp, sina, sinp, csca, cscp, cotp, cota, nbpa, npam, cosb, C1, C2, F1, F2, F3, G1, G2, G3, G4, H1, H2, H3, H4, GradAcosb] = outer_dihedral_hessian_p_terms();

            coord3d GradAcosa = (aph * cosa - amh) * rapn + (amh * cosa - aph) * ramn;
            coord3d GradAsina = -cosa * csca * GradAcosa;
            coord3d GradAcsca = -GradAsina * csca * csca;
            coord3d GradAcota = (sina * GradAcosa - cosa * GradAsina) * csca * csca;
            coord3d GradAcosp = (pbh - pah * cosp) * rapn;
            coord3d GradAsinp = -cosp * cscp * GradAcosp;
            coord3d GradAcscp = -GradAsinp * cscp * cscp;
            coord3d GradAcotp = (sinp * GradAcosp - cosp * GradAsinp) * cscp * cscp;

            mat3 GradAamh = (tensor_product(amh, amh) - identity3()) * ramn;
            coord3d GradAramn = amh * ramn * ramn;
            mat3 GradAaph = (tensor_product(aph, aph) - identity3()) * rapn;
            coord3d GradArapn = aph * rapn * rapn;

            mat3 GradApah = (identity3() - tensor_product(pah, pah)) * rapn;
            coord3d GradArpan = -pah * rapn * rapn;

            coord3d GradAC1 = cscp * GradArpan + rapn * GradAcscp;
            coord3d GradAC2 = cosb * (cota * GradAcsca + csca * GradAcota) + cota * csca * GradAcosb;

            mat3 GradAnpam = (sina * (cross(aph, GradAamh) - cross(amh, GradAaph)) - tensor_product(cross(aph, amh), GradAsina)) * csca * csca;

            mat3 GradAF1 = -cross(pbh, GradAnpam);
            mat3 GradAF2 = tensor_product(pah, sinp * GradAcosb + cosb * GradAsinp) + GradApah * sinp * cosb;
            mat3 GradAF3 = cotp * cosb * (-tensor_product(pah, GradAcosp) - GradApah * cosp) + tensor_product(pbh - pah * cosp, GradAcotp * cosb) + tensor_product(pbh - pah * cosp, cotp * GradAcosb);

            mat3 GradAnbpa = (sinp * (cross(pbh, GradApah)) - tensor_product(cross(pbh, pah), GradAsinp)) * cscp * cscp;

            mat3 GradAG1 = tensor_product(aph, GradAcosb) + GradAaph * cosb;
            mat3 GradAG2 = tensor_product(cross(amh, nbpa), GradAcsca) + csca * (cross(amh, GradAnbpa) - cross(nbpa, GradAamh));
            mat3 GradAG3 = tensor_product(amh, GradAcosb) + GradAamh * cosb;
            mat3 GradAG4 = tensor_product(cross(nbpa, aph), GradAcsca) + csca * (cross(nbpa, GradAaph) - cross(aph, GradAnbpa));

            mat3 GradAH1 = tensor_product(aph, GradAcosa * rapn + GradArapn * cosa) + GradAaph * cosa * rapn;
            mat3 GradAH2 = tensor_product(amh, GradArapn) + GradAamh * rapn;
            mat3 GradAH3 = tensor_product(amh, GradAcosa * ramn + GradAramn * cosa) + GradAamh * cosa * ramn;
            mat3 GradAH4 = tensor_product(aph, GradAramn) + GradAaph * ramn;

            mat3 GradGradAcosb = C1 * (GradAF1 - GradAF2 + GradAF3) + tensor_product(F1 - F2 + F3, GradAC1) + rapn * (GradAG1 - GradAG2) + tensor_product(G1 - G2, GradArapn) + ramn * (GradAG3 - GradAG4) + tensor_product(G3 - G4, GradAramn) + C2 * (GradAH1 - GradAH2 + GradAH3 - GradAH4) + tensor_product(H1 - H2 + H3 - H4, GradAC2);
            return d_get(c.f_outer_dihedral, j) * harmonic_energy_hessian(d_get(c.outer_dih0_p, j), cosb, GradAcosb, GradAcosb, GradGradAcosb); // Harmonic Energy Hessian
        }

        mat3 outer_dihedral_hessian_p_b(const Constants &c) const
        {
            auto [cosa, cosp, sina, sinp, csca, cscp, cotp, cota, nbpa, npam, cosb, C1, C2, F1, F2, F3, G1, G2, G3, G4, H1, H2, H3, H4, GradAcosb] = outer_dihedral_hessian_p_terms();

            coord3d GradBcosp = (pah - pbh * cosp) * rbpn;
            coord3d GradBsinp = -cosp * cscp * GradBcosp;
            coord3d GradBcscp = -GradBsinp * cscp * cscp;
            coord3d GradBcotp = (sinp * GradBcosp - cosp * GradBsinp) * cscp * cscp;

            mat3 GradBpbh = (identity3() - tensor_product(pbh, pbh)) * rbpn;
            coord3d GradBrpbn = -pbh * rbpn * rbpn;
            coord3d GradBcosb = ((dot(npam, cross(pah, pbh)) * pbh - cross(npam, pah)) * sinp * rbpn - dot(npam, cross(pbh, pah)) * GradBsinp) * cscp * cscp;

            coord3d GradBC1 = GradBcscp * rapn;
            coord3d GradBC2 = cota * csca * GradBcosb;

            mat3 GradBF1 = cross(npam, GradBpbh);
            mat3 GradBF2 = tensor_product(pah, sinp * GradBcosb + cosb * GradBsinp);
            mat3 GradBF3 = cotp * cosb * (GradBpbh - tensor_product(pah, GradBcosp)) + tensor_product(pbh - pah * cosp, GradBcotp * cosb) + tensor_product(pbh - pah * cosp, cotp * GradBcosb);

            mat3 GradBnbpa = (sinp * (-cross(pah, GradBpbh)) - tensor_product(cross(pbh, pah), GradBsinp)) * cscp * cscp;
            mat3 GradBG1 = tensor_product(aph, GradBcosb);
            mat3 GradBG2 = csca * (cross(amh, GradBnbpa));
            mat3 GradBG3 = tensor_product(amh, GradBcosb);
            mat3 GradBG4 = csca * (-cross(aph, GradBnbpa));

            mat3 GradGradBcosb = C1 * (GradBF1 - GradBF2 + GradBF3) + tensor_product(F1 - F2 + F3, GradBC1) + rapn * (GradBG1 - GradBG2) + ramn * (GradBG3 - GradBG4) + tensor_product(H1 - H2 + H3 - H4, GradBC2);

            return d_get(c.f_outer_dihedral, j) * harmonic_energy_hessian(d_get(c.outer_dih0_p, j), cosb, GradAcosb, GradBcosb, GradGradBcosb); // Harmonic Energy Hessian
        }

        mat3 outer_dihedral_hessian_p_m(const Constants &c) const
        {
            auto [cosa, cosp, sina, sinp, csca, cscp, cotp, cota, nbpa, npam, cosb, C1, C2, F1, F2, F3, G1, G2, G3, G4, H1, H2, H3, H4, GradAcosb] = outer_dihedral_hessian_p_terms();

            coord3d GradMcosa = (aph - amh * cosa) * ramn;
            coord3d GradMsina = -cosa * csca * GradMcosa;
            coord3d GradMcsca = -GradMsina * csca * csca;
            coord3d GradMcota = (sina * GradMcosa - cosa * GradMsina) * csca * csca;
            coord3d GradMcosb = ((cross(nbpa, aph) - dot(cross(nbpa, aph), amh) * amh) * sina * ramn - dot(nbpa, cross(aph, amh)) * GradMsina) * csca * csca;

            mat3 GradMamh = (identity3() - tensor_product(amh, amh)) * ramn;
            coord3d GradMramn = -amh * ramn * ramn;

            coord3d GradMC2 = cosb * (cota * GradMcsca + csca * GradMcota) + cota * csca * GradMcosb;
            mat3 GradMnpam = (sina * (cross(aph, GradMamh)) - tensor_product(cross(aph, amh), GradMsina)) * csca * csca;
            mat3 GradMF1 = -cross(pbh, GradMnpam);
            mat3 GradMF2 = tensor_product(pah, sinp * GradMcosb);
            mat3 GradMF3 = tensor_product(pbh - pah * cosp, cotp * GradMcosb);

            mat3 GradMG1 = tensor_product(aph, GradMcosb);
            mat3 GradMG2 = tensor_product(cross(amh, nbpa), GradMcsca) + csca * (-cross(nbpa, GradMamh));
            mat3 GradMG3 = tensor_product(amh, GradMcosb) + GradMamh * cosb;
            mat3 GradMG4 = tensor_product(cross(nbpa, aph), GradMcsca);

            mat3 GradMH1 = tensor_product(aph, GradMcosa * rapn);
            mat3 GradMH2 = GradMamh * rapn;
            mat3 GradMH3 = tensor_product(amh, GradMcosa * ramn + GradMramn * cosa) + GradMamh * cosa * ramn;
            mat3 GradMH4 = tensor_product(aph, GradMramn);

            mat3 GradGradMcosb = C1 * (GradMF1 - GradMF2 + GradMF3) + rapn * (GradMG1 - GradMG2) + ramn * (GradMG3 - GradMG4) + tensor_product(G3 - G4, GradMramn) + C2 * (GradMH1 - GradMH2 + GradMH3 - GradMH4) + tensor_product(H1 - H2 + H3 - H4, GradMC2);
            return d_get(c.f_outer_dihedral, j) * harmonic_energy_hessian(d_get(c.outer_dih0_p, j), cosb, GradAcosb, GradMcosb, GradGradMcosb); // Harmonic Energy Hessian
        }

        mat3 outer_dihedral_hessian_p_p(const Constants &c) const
        {
            auto [cosa, cosp, sina, sinp, csca, cscp, cotp, cota, nbpa, npam, cosb, C1, C2, F1, F2, F3, G1, G2, G3, G4, H1, H2, H3, H4, GradAcosb] = outer_dihedral_hessian_p_terms();

            coord3d GradPcosa = (amh - aph * cosa) * rapn;
            coord3d GradPsina = -cosa * csca * GradPcosa;
            coord3d GradPcsca = -GradPsina * csca * csca;
            coord3d GradPcota = (sina * GradPcosa - cosa * GradPsina) * csca * csca;
            coord3d GradPcosp = (pbh * cosp - pah) * rbpn + (pah * cosp - pbh) * rapn;
            coord3d GradPsinp = -cosp * cscp * GradPcosp;
            coord3d GradPcscp = -GradPsinp * cscp * cscp;
            coord3d GradPcotp = (sinp * GradPcosp - cosp * GradPsinp) * cscp * cscp;

            mat3 GradPaph = (identity3() - tensor_product(aph, aph)) * rapn;
            coord3d GradPrapn = -aph * rapn * rapn;
            mat3 GradPpbh = (tensor_product(pbh, pbh) - identity3()) * rbpn;
            coord3d GradPrpbn = pbh * rbpn * rbpn;
            mat3 GradPpah = (tensor_product(pah, pah) - identity3()) * rapn;
            coord3d GradPrpan = pah * rapn * rapn;
            coord3d GradPcosbP1 = (((dot(cross(npam, pbh), pah) * pah - cross(npam, pbh)) * rapn - (dot(cross(npam, pah), pbh) * pbh - cross(npam, pah)) * rbpn) * sinp - dot(npam, cross(pbh, pah)) * GradPsinp) * cscp * cscp;
            coord3d GradPcosbP2 = (-(cross(nbpa, amh) - dot(cross(nbpa, amh), aph) * aph) * sina * rapn - dot(nbpa, cross(aph, amh)) * GradPsina) * csca * csca;
            coord3d GradPcosb = GradPcosbP1 + GradPcosbP2;

            coord3d GradPC1 = cscp * GradPrpan + rapn * GradPcscp;
            coord3d GradPC2 = cosb * (cota * GradPcsca + csca * GradPcota) + cota * csca * GradPcosb;

            mat3 GradPnpam = (sina * (-cross(amh, GradPaph)) - tensor_product(cross(aph, amh), GradPsina)) * csca * csca;
            mat3 GradPF1 = cross(npam, GradPpbh) - cross(pbh, GradPnpam);
            mat3 GradPF2 = tensor_product(pah, sinp * GradPcosb + cosb * GradPsinp) + GradPpah * sinp * cosb;
            mat3 GradPF3 = cotp * cosb * (GradPpbh - tensor_product(pah, GradPcosp) - GradPpah * cosp) + tensor_product(pbh - pah * cosp, GradPcotp * cosb) + tensor_product(pbh - pah * cosp, cotp * GradPcosb);

            mat3 GradPG1 = tensor_product(aph, GradPcosb) + GradPaph * cosb;
            mat3 GradPnbpa = (sinp * (cross(pbh, GradPpah) - cross(pah, GradPpbh)) - tensor_product(cross(pbh, pah), GradPsinp)) * cscp * cscp;
            mat3 GradPG2 = tensor_product(cross(amh, nbpa), GradPcsca) + csca * (cross(amh, GradPnbpa));
            mat3 GradPG3 = tensor_product(amh, GradPcosb);
            mat3 GradPG4 = tensor_product(cross(nbpa, aph), GradPcsca) + csca * (cross(nbpa, GradPaph) - cross(aph, GradPnbpa));

            mat3 GradPH1 = tensor_product(aph, GradPcosa * rapn + GradPrapn * cosa) + GradPaph * cosa * rapn;
            mat3 GradPH2 = tensor_product(amh, GradPrapn);
            mat3 GradPH3 = tensor_product(amh, GradPcosa * ramn);
            mat3 GradPH4 = GradPaph * ramn;

            mat3 GradGradPcosb = C1 * (GradPF1 - GradPF2 + GradPF3) + tensor_product(F1 - F2 + F3, GradPC1) + rapn * (GradPG1 - GradPG2) + tensor_product(G1 - G2, GradPrapn) + ramn * (GradPG3 - GradPG4) + C2 * (GradPH1 - GradPH2 + GradPH3 - GradPH4) + tensor_product(H1 - H2 + H3 - H4, GradPC2);
            return d_get(c.f_outer_dihedral, j) * harmonic_energy_hessian(d_get(c.outer_dih0_p, j), cosb, GradAcosb, GradPcosb, GradGradPcosb); // Harmonic Energy Hessian
        }

        // Computes gradient related to bending of outer angles. ~20 FLOPs
        /**
         * @brief Compute the gradient of the outer angle-m term.
         * @param c The constants for the threadIdx^th node.
         * @return The gradient of the outer angle-m term.
         */
        coord3d outer_angle_gradient_m(const Constants &c) const
        {
            real_t cos_angle = -dot(abh, bmh);                                                                          // Compute outer angle. ab,bm
            coord3d grad = (bmh + abh * cos_angle) * rabn;                                                              // Derivative of outer angles Eq. 30. Buster Thesis
            return d_get(c.f_outer_angle_m, j) * harmonic_energy_gradient(d_get(c.outer_angle_m0, j), cos_angle, grad); // Harmonic Energy Gradient: Eq. 30 multiplied by harmonic term.
        }

        /**
         * @brief Compute the gradient of the outer angle-p term.
         * @param c The constants for the threadIdx^th node.
         * @return The gradient of the outer angle-p term.
         */
        coord3d outer_angle_gradient_p(const Constants &c) const
        {
            real_t cos_angle = -dot(abh, bph);                                                                          // Compute outer angle. ab,bp
            coord3d grad = (bph + abh * cos_angle) * rabn;                                                              // Derivative of outer angles Eq. 28. Buster Thesis
            return d_get(c.f_outer_angle_p, j) * harmonic_energy_gradient(d_get(c.outer_angle_p0, j), cos_angle, grad); // Harmonic Energy Gradient: Eq. 28 multiplied by harmonic term.
        }
        // Chain rule terms for dihedral calculation
        // Computes gradient related to dihedral/out-of-plane term. ~75 FLOPs
        /**
         * @brief Compute the gradient of the inner dihedral term.
         * @param c The constants for the threadIdx^th node.
         * @return The gradient of the inner dihedral term.
         */
        coord3d inner_dihedral_gradient(const Constants &c) const
        {
            coord3d nabc, nbcd;
            real_t cos_b, cos_c, r_sin_b, r_sin_c;
            cos_b = dot(bah, bch);
            r_sin_b = sycl::rsqrt((real_t)1.0 - cos_b * cos_b);
            nabc = cross(bah, bch) * r_sin_b;
            cos_c = dot(-bch, cdh);
            r_sin_c = sycl::rsqrt((real_t)1.0 - cos_c * cos_c);
            nbcd = cross(-bch, cdh) * r_sin_c;

            real_t cos_beta = dot(nabc, nbcd);        // Inner dihedral angle from planes abc,bcd.
            real_t cot_b = cos_b * r_sin_b * r_sin_b; // cos(b)/sin(b)^2

            // Derivative w.r.t. inner dihedral angle F and G in Eq. 26
            coord3d grad = cross(bch, nbcd) * r_sin_b * rabn - bah * cos_beta * rabn + (cot_b * cos_beta * rabn) * (bch - bah * cos_b);
            return d_get(c.f_inner_dihedral, j) * harmonic_energy_gradient(d_get(c.inner_dih0, j), cos_beta, grad); // Eq. 26.
        }

        // Computes gradient from dihedral angles constituted by the planes bam, amp ~162 FLOPs
        /**
         * @brief Compute the gradient of the outer dihedral-a term.
         * @param c The constants for the threadIdx^th node.
         * @return The gradient of the outer dihedral-a term.
         */
        coord3d outer_dihedral_gradient_a(const Constants &c) const
        {
            coord3d nbam_hat, namp_hat;
            real_t cos_a, cos_m, r_sin_a, r_sin_m;

            cos_a = dot(abh, amh);
            r_sin_a = sycl::rsqrt((real_t)1.0 - cos_a * cos_a);
            nbam_hat = cross(abh, amh) * r_sin_a;
            cos_m = dot(-amh, mph);
            r_sin_m = sycl::rsqrt((real_t)1.0 - cos_m * cos_m);
            namp_hat = cross(-amh, mph) * r_sin_m;

            real_t cos_beta = dot(nbam_hat, namp_hat); // Outer Dihedral angle bam, amp
            real_t cot_a = cos_a * r_sin_a * r_sin_a;
            real_t cot_m = cos_m * r_sin_m * r_sin_m;

            // Derivative w.r.t. outer dihedral angle, factorized version of Eq. 31.
            coord3d grad = cross(mph, nbam_hat) * ramn * r_sin_m - (cross(namp_hat, abh) * ramn + cross(amh, namp_hat) * rabn) * r_sin_a +
                           cos_beta * (abh * rabn + ramn * ((real_t)2.0 * amh + cot_m * (mph + cos_m * amh)) - cot_a * (ramn * (abh - amh * cos_a) + rabn * (amh - abh * cos_a)));

            // Eq. 31 multiplied by harmonic term.
            return d_get(c.f_outer_dihedral, j) * harmonic_energy_gradient(d_get(c.outer_dih0_a, j), cos_beta, grad);
        }

        // Computes gradient from dihedral angles constituted by the planes nbmp, nmpa ~92 FLOPs
        /**
         * @brief Compute the gradient of the outer dihedral-m term.
         * @param c The constants for the threadIdx^th node.
         * @return The gradient of the outer dihedral-m term.
         */
        coord3d outer_dihedral_gradient_m(const Constants &c) const
        {
            coord3d nbmp_hat, nmpa_hat;
            real_t cos_m, cos_p, r_sin_m, r_sin_p;
            cos_m = dot(mbh, mph);
            r_sin_m = sycl::rsqrt((real_t)1.0 - cos_m * cos_m);
            nbmp_hat = cross(mbh, mph) * r_sin_m;
            cos_p = dot(-mph, pah);
            r_sin_p = sycl::rsqrt((real_t)1.0 - cos_p * cos_p);
            nmpa_hat = cross(-mph, pah) * r_sin_p;

            // Cosine to the outer dihedral angle constituted by the planes bmp and mpa
            real_t cos_beta = dot(nbmp_hat, nmpa_hat); // Outer dihedral angle bmp,mpa.
            real_t cot_p = cos_p * r_sin_p * r_sin_p;

            // Derivative w.r.t. outer dihedral angle, factorized version of Eq. 32.
            coord3d grad = rapn * (cot_p * cos_beta * (-mph - pah * cos_p) - cross(nbmp_hat, mph) * r_sin_p - pah * cos_beta);

            // Eq. 32 multiplied by harmonic term.
            return d_get(c.f_outer_dihedral, j) * harmonic_energy_gradient(d_get(c.outer_dih0_m, j), cos_beta, grad);
        }

        // Computes gradient from dihedral angles constituted by the planes bpa, pam ~162 FLOPs
        /**
         * @brief Compute the gradient of the outer dihedral-p term.
         * @param c The constants for the threadIdx^th node.
         * @return The gradient of the outer dihedral-p term.
         */
        coord3d outer_dihedral_gradient_p(const Constants &c) const
        {
            coord3d nbpa_hat, npam_hat;
            real_t cos_p, cos_a, r_sin_p, r_sin_a;
            cos_a = dot(aph, amh);
            r_sin_a = sycl::rsqrt((real_t)1.0 - cos_a * cos_a);
            npam_hat = cross(aph, amh) * r_sin_a;
            cos_p = dot(pbh, -aph);
            r_sin_p = sycl::rsqrt((real_t)1.0 - cos_p * cos_p);
            nbpa_hat = cross(pbh, -aph) * r_sin_p;

            real_t cos_beta = dot(nbpa_hat, npam_hat); // Outer dihedral angle bpa, pam.
            real_t cot_p = cos_p * r_sin_p * r_sin_p;
            real_t cot_a = cos_a * r_sin_a * r_sin_a;

            // Derivative w.r.t. outer dihedral angle, factorized version of Eq. 33.
            coord3d grad = cross(npam_hat, pbh) * rapn * r_sin_p - (cross(amh, nbpa_hat) * rapn + cross(nbpa_hat, aph) * ramn) * r_sin_a +
                           cos_beta * (amh * ramn + rapn * ((real_t)2.0 * aph + cot_p * (pbh + cos_p * aph)) - cot_a * (rapn * (amh - aph * cos_a) + ramn * (aph - amh * cos_a)));

            // Eq. 33 multiplied by harmonic term.
            return d_get(c.f_outer_dihedral, j) * harmonic_energy_gradient(d_get(c.outer_dih0_p, j), cos_beta, grad);
        }

        // Internal coordinate gradients
        /**
         * @brief Compute the gradient of the bond length term.
         * @param c The constants for the threadIdx^th node.
         * @return The gradient of the bond length term.
         */
        coord3d bond_length_gradient(const Constants &c) const
        {
            return c.f_bond[j] * harmonic_energy_gradient(bond(), c.r0[j], abh);
        }
        // Sum of angular gradient components.
        /**
         * @brief Compute the sum of the gradients of the bending terms.
         * @param c The constants for the threadIdx^th node.
         * @return The sum of the gradients of the bending terms.
         */
        coord3d angle_gradient(const Constants &c) const { return inner_angle_gradient(c) + outer_angle_gradient_p(c) + outer_angle_gradient_m(c); }
        // Sum of inner and outer dihedral gradient components.
        /**
         * @brief Compute the sum of the gradients of the dihedral terms.
         * @param c The constants for the threadIdx^th node.
         * @return The sum of the gradients of the dihedral terms.
         */
        coord3d dihedral_gradient(const Constants &c) const
        {
            switch (FFT)
            {
            case PEDERSEN:
                return inner_dihedral_gradient(c) + outer_dihedral_gradient_a(c) + outer_dihedral_gradient_m(c) + outer_dihedral_gradient_p(c);
            case WIRZ:
                return inner_dihedral_gradient(c);
            default:
                return inner_dihedral_gradient(c) + outer_dihedral_gradient_a(c) + outer_dihedral_gradient_m(c) + outer_dihedral_gradient_p(c);
            }
        }
        /**
         * @brief Compute the energy contribution of the bond length term.
         * @param c The constants for the threadIdx^th node.
         * @return The energy contribution of the bond length term.
         */
        real_t bond_energy(const Constants &c) const
        {
            return (real_t)0.5 * d_get(c.f_bond, j) * harmonic_energy(bond(), d_get(c.r0, j));
        }
        /**
         * @brief Compute the total energy contribution of the bending terms.
         * @param c The constants for the threadIdx^th node.
         * @return The energy contribution of the bending terms.
         */
        real_t bend_energy(const Constants &c) const
        {
            return d_get(c.f_inner_angle, j) * harmonic_energy(angle(), d_get(c.angle0, j));
        }

        /**
         * @brief Compute the total energy contribution of the dihedral terms.
         * @param c The constants for the threadIdx^th node.
         * @return The energy contribution of the dihedral terms.
         */
        real_t dihedral_energy(const Constants &c) const
        {
            return d_get(c.f_inner_dihedral, j) * harmonic_energy(dihedral(), d_get(c.inner_dih0, j));
        }
        // Harmonic energy contribution from bond stretching, angular bending and dihedral angle bending.
        // 71 FLOPs
        /**
         * @brief Compute the total energy contribution of the bond length, bending and dihedral terms.
         * @param c The constants for the threadIdx^th node.
         * @return The energy contribution of the bond length, bending and dihedral terms.
         */
        real_t energy(const Constants &c) const
        {
            switch (FFT)
            {
            case FLAT_BOND:
                return bond_energy(c);
            default:
                return bond_energy(c) + bend_energy(c) + dihedral_energy(c);
            }
        }
        // Sum of bond, angular and dihedral gradient components.
        /**
         * @brief Compute the total gradient of the bond length, bending and dihedral terms w.r.t. the coordinates of the threadIdx^th node.
         * @param c The constants for the threadIdx^th node.
         * @return The gradient of the bond length, bending and dihedral terms w.r.t. the coordinates of the threadIdx^th node.
         */
        coord3d gradient(const Constants &c) const
        {
            switch (FFT)
            {
            case FLAT_BOND:
                return bond_length_gradient(c);
            case BOND:
                return bond_length_gradient(c);
            case ANGLE:
                return inner_angle_gradient(c);
            case DIH:
                return inner_dihedral_gradient(c);
            case ANGLE_M:
                return outer_angle_gradient_m(c);
            case ANGLE_P:
                return outer_angle_gradient_p(c);
            case DIH_A:
                return outer_dihedral_gradient_a(c);
            case DIH_M:
                return outer_dihedral_gradient_m(c);
            case DIH_P:
                return outer_dihedral_gradient_p(c);
            default:
                return bond_length_gradient(c) + angle_gradient(c) + dihedral_gradient(c);
            }
        }

        mat3 hessian_a(const Constants &c) const
        {
            switch (FFT)
            {
            case BOND:
                return bond_hessian_a(c);
            case ANGLE:
                return inner_angle_hessian_a(c);
            case ANGLE_M:
                return outer_angle_hessian_m_a(c);
            case ANGLE_P:
                return outer_angle_hessian_p_a(c);
            case DIH:
                return dihedral_hessian_a(c);
            case DIH_A:
                return outer_dihedral_hessian_a_a(c);
            case DIH_M:
                return outer_dihedral_hessian_m_a(c);
            case DIH_P:
                return outer_dihedral_hessian_p_a(c);
            default:
                return bond_hessian_a(c) + inner_angle_hessian_a(c) + outer_angle_hessian_m_a(c) + outer_angle_hessian_p_a(c) + dihedral_hessian_a(c) + outer_dihedral_hessian_a_a(c) + outer_dihedral_hessian_m_a(c) + outer_dihedral_hessian_p_a(c);
            }
        }

        mat3 hessian_b(const Constants &c) const
        {
            switch (FFT)
            {
            case BOND:
                return bond_hessian_b(c);
            case ANGLE:
                return inner_angle_hessian_b(c);
            case ANGLE_M:
                return outer_angle_hessian_m_b(c);
            case ANGLE_P:
                return outer_angle_hessian_p_b(c);
            case DIH:
                return dihedral_hessian_b(c);
            case DIH_A:
                return outer_dihedral_hessian_a_b(c);
            case DIH_M:
                return outer_dihedral_hessian_m_b(c);
            case DIH_P:
                return outer_dihedral_hessian_p_b(c);
            default:
                return bond_hessian_b(c) + inner_angle_hessian_b(c) + outer_angle_hessian_m_b(c) + outer_angle_hessian_p_b(c) + dihedral_hessian_b(c) + outer_dihedral_hessian_a_b(c) + outer_dihedral_hessian_m_b(c) + outer_dihedral_hessian_p_b(c);
            }
        }

        mat3 hessian_c(const Constants &c) const
        {
            switch (FFT)
            {
            case BOND:
                return mat3();
            case ANGLE:
                return inner_angle_hessian_c(c);
            case ANGLE_M:
                return mat3();
            case ANGLE_P:
                return mat3();
            case DIH:
                return dihedral_hessian_c(c);
            case DIH_A:
                return mat3();
            case DIH_M:
                return mat3();
            case DIH_P:
                return mat3();
            default:
                return inner_angle_hessian_c(c) + dihedral_hessian_c(c);
            }
        }

        mat3 hessian_d(const Constants &c) const
        {
            switch (FFT)
            {
            case BOND:
                return mat3();
            case ANGLE:
                return mat3();
            case ANGLE_M:
                return mat3();
            case ANGLE_P:
                return mat3();
            case DIH:
                return dihedral_hessian_d(c);
            case DIH_A:
                return mat3();
            case DIH_M:
                return mat3();
            case DIH_P:
                return mat3();
            default:
                return dihedral_hessian_d(c);
            }
            return dihedral_hessian_d(c);
        }

        mat3 hessian_m(const Constants &c) const
        {
            switch (FFT)
            {
            case BOND:
                return mat3();
            case ANGLE:
                return mat3();
            case ANGLE_M:
                return outer_angle_hessian_m_m(c);
            case ANGLE_P:
                return mat3();
            case DIH:
                return mat3();
            case DIH_A:
                return outer_dihedral_hessian_a_m(c);
            case DIH_M:
                return outer_dihedral_hessian_m_m(c);
            case DIH_P:
                return outer_dihedral_hessian_p_m(c);
            default:
                return outer_angle_hessian_m_m(c) + outer_dihedral_hessian_a_m(c) + outer_dihedral_hessian_m_m(c) + outer_dihedral_hessian_p_m(c);
            }
        }

        mat3 hessian_p(const Constants &c) const
        {
            switch (FFT)
            {
            case BOND:
                return mat3();
            case ANGLE:
                return mat3();
            case ANGLE_M:
                return mat3();
            case ANGLE_P:
                return outer_angle_hessian_p_p(c);
            case DIH:
                return mat3();
            case DIH_A:
                return outer_dihedral_hessian_a_p(c);
            case DIH_M:
                return outer_dihedral_hessian_m_p(c);
            case DIH_P:
                return outer_dihedral_hessian_p_p(c);
            default:
                return outer_angle_hessian_p_p(c) + outer_dihedral_hessian_a_p(c) + outer_dihedral_hessian_m_p(c) + outer_dihedral_hessian_p_p(c);
            }
        }

        // Reciprocal lengths of arcs ab, ac, am, ap.
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

        // Base Arcs,
        coord3d
            ab,
            ac,
            ad;

        // All normalized arcs required to perform energy & gradient calculations.
        // Note that all these arcs are cyclical the arc ab becomes: ab->ac->ad,  the arc ac becomes: ac->ad->ab , the arc bc becomes: bc->cd->db (For iterations 0, 1, 2)
        // As such the naming convention here is related to the arcs as they are used in the 0th iteration.
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
            face_center, // Center of the face to the left of the arc a->b, a->b, a->c
            face_offset; // Difference between the node coordinates X and the face-center coordinates face_center

        coord3d A[3];
    };

    /**
     * @brief Compute the total gradient of the bond, flatness, bending and dihedral terms w.r.t. the coordinates of the threadIdx^th node.
     * @param c The constants for the threadIdx^th node.
     * @return The gradient of the bond, flatness, bending and dihedral terms w.r.t. the coordinates of the threadIdx^th node.
     */
    coord3d gradient(const Span<coord3d> X) const
    {
        sycl::group_barrier(cta);
        coord3d grad = {0.0, 0.0, 0.0};
        for (int j = 0; j < 3; j++)
        {
            ArcData arc = ArcData(cta.get_local_linear_id(), j, X, node_graph);
            grad += arc.gradient(constants);
        }
        switch (FFT)
        {
         case FLATNESS_ENABLED: {
             FaceData face(cta, X, node_graph);
             auto face_grad = face.flatness_gradient(constants, reinterpret_cast<coord3d*>(sdata));
             return grad + face_grad;
             }
         case FLAT_BOND:{
             FaceData face(cta, X, node_graph);
             auto face_grad = face.flatness_gradient(constants, reinterpret_cast<coord3d*>(sdata));
             return grad + face_grad;
             }
         default:{
             return grad;
            }
         }
    }

    hessian_t<T, K> hessian(const Span<coord3d> X) const
    {
        sycl::group_barrier(cta);
        hessian_t<T, K> hess(node_graph, node_id);
        for (int j = 0; j < 3; j++)
        {
            ArcData arc = ArcData(cta.get_local_linear_id(), j, X, node_graph);
            hess.A[0] += arc.hessian_a(constants);
            hess.A[1 + j] += arc.hessian_b(constants);
            hess.A[1 + (j + 1) % 3] += arc.hessian_c(constants);
            hess.A[1 + (j + 2) % 3] += arc.hessian_d(constants);
            hess.A[4 + j] += arc.hessian_m(constants);
            hess.A[7 + j] += arc.hessian_p(constants);
        }
        return hess;
    }

    // Uses finite difference to compute the hessian
    hessian_t<T, K> fd_hessian(const Span<coord3d> X, const float reldelta = 1e-7) const
    {
        hessian_t<T, K> hess_fd(node_graph, node_id);
        for (uint16_t i = 0; i < N; i++)
        {
            for (int j = 0; j < 10; j++)
            {
                auto node = hess_fd.indices[j];
                coord3d X0 = X[node];
                for (int k = 0; k < 3; k++)
                {
                    if (i == node_id)
                    {
                        X[node][k] = X0[k] + X0[k] * reldelta;
                    }
                    coord3d grad_X0_p = gradient(X);
                    sycl::group_barrier(cta);
                    if (i == node_id)
                    {
                        X[node][k] = X0[k] - X0[k] * reldelta;
                    }
                    coord3d grad_X0_m = gradient(X);
                    sycl::group_barrier(cta);
                    if (i == node_id)
                    {
                        hess_fd.A[j][0][k] = (grad_X0_p[0] - grad_X0_m[0]) / (2 * X0[k] * reldelta);
                        hess_fd.A[j][1][k] = (grad_X0_p[1] - grad_X0_m[1]) / (2 * X0[k] * reldelta);
                        hess_fd.A[j][2][k] = (grad_X0_p[2] - grad_X0_m[2]) / (2 * X0[k] * reldelta);
                        X[node][k] = X0[k];
                    }
                    sycl::group_barrier(cta);
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
    real_t energy(const Span<coord3d> X) const
    {
        sycl::group_barrier(cta);
        real_t arc_energy = (real_t)0.0;

        //(71 + 124) * 3 * N  = 585*N FLOPs
        for (uint8_t j = 0; j < 3; j++)
        {
            ArcData arc = ArcData(cta.get_local_linear_id(), j, X, node_graph);
            arc_energy += arc.energy(constants);
        }
        return sycl::reduce_over_group(cta, arc_energy, sycl::plus<real_t>{});

        // switch (FFT)
        //{
        // case FLATNESS_ENABLED: {
        //     FaceData face(cta, X, node_graph);
        //     return reduce_over_group(cta, arc_energy + face.flatness_energy(constants), sycl::plus<real_t>{});
        //     }
        // case FLAT_BOND: {
        //     FaceData face(cta, X, node_graph);
        //     return reduce_over_group(cta, arc_energy + face.flatness_energy(constants), sycl::plus<real_t>{});
        //     }
        // default:
        //     return reduce_over_group(cta, arc_energy, sycl::plus<real_t>{});
        // }
    }

    // Golden Section Search, using fixed iterations.
    /**
     * @brief Golden Section Search for line-search.
     * @param X The coordinates of the nodes.
     * @param r0 The direction of the line-search.
     * @param X1 memory for storing temporary coordinates at x1.
     * @param X2 memory for storing temporary coordinates at x2.
     * @return The step-size alpha
     */
    real_t GSS(const Span<coord3d> X, const coord3d &r0, const Span<coord3d> X1, const Span<coord3d> X2) const
    {
        const real_t tau = (real_t)0.6180339887;
        // Line search x - values;
        real_t a = 0.0;
        real_t b = (real_t)1.0;

        real_t x1, x2;
        x1 = (a + ((real_t)1. - tau) * (b - a));
        x2 = (a + tau * (b - a));
        // Actual coordinates resulting from each traversal
        X1[node_id] = X[node_id] + x1 * r0;
        X2[node_id] = X[node_id] + x2 * r0;

        real_t f1 = energy(X1);
        real_t f2 = energy(X2);

        for (int i = 0; i < 20; i++)
        {
            if (f1 > f2)
            {
                a = x1;
                x1 = x2;
                f1 = f2;
                x2 = a + tau * (b - a);
                X2[node_id] = X[node_id] + x2 * r0;
                f2 = energy(X2);
            }
            else
            {
                b = x2;
                x2 = x1;
                f2 = f1;
                x1 = a + ((real_t)1.0 - tau) * (b - a);
                X1[node_id] = X[node_id] + x1 * r0;
                f1 = energy(X1);
            }
        }
        if (f1 > energy(X))
        {
            return (real_t)0.0;
        }
        // Line search coefficient
        real_t alpha = (a + b) / (real_t)2.0;
        return alpha;
    }

    /**
     * @brief Conjugate Gradient Method for energy minimization.
     * @param X The coordinates of the nodes.
     * @param X1 memory for storing temporary coordinates.
     * @param X2 memory for storing temporary coordinates.
     * @param MaxIter The maximum number of iterations.
     */
    void CG(const Span<coord3d> X, const Span<coord3d> X1, const Span<coord3d> X2, const size_t MaxIter)
    {
        real_t alpha, beta, g0_norm2, s_norm;
        coord3d g0, g1, s;
        g0 = gradient(X);
        s = -g0;

        // Normalize To match reference python implementation by Buster.
        s_norm = SQRT(sycl::reduce_over_group(cta, dot(s, s), sycl::plus<real_t>{}));
        // s_norm = SQRT(reduction(sdata, dot(s,s)));
        s /= s_norm;

        sycl::group_barrier(cta);
        for (size_t i = 0; i < MaxIter; i++)
        {
            alpha = GSS(X, s, X1, X2);

            if (alpha > (real_t)0.0)
            {
                X1[node_id] = X[node_id] + alpha * s;
            }
            g1 = gradient(X1);

            // Polak Ribiere method
            g0_norm2 = sycl::reduce_over_group(cta, dot(g0, g0), sycl::plus<real_t>{});
            beta = sycl::max(sycl::reduce_over_group(cta, dot(g1, (g1 - g0)), sycl::plus<real_t>{}) / g0_norm2, (real_t)0.0);

            if (alpha > (real_t)0.0)
            {
                X[node_id] = X1[node_id];
            }
            else
            {
                g1 = g0;
                beta = (real_t)0.0;
            }

            s = -g1 + beta * s;
            g0 = g1;
// Normalize Search Direction using MaxNorm or 2Norm
#if USE_MAX_NORM == 1
            s_norm = sycl::reduce_over_group(cta, sycl::max(sycl::max(s.x, s.y), s.z), sycl::greater<real_t>{});
#else
            s_norm = SQRT(sycl::reduce_over_group(cta, dot(s, s), sycl::plus<real_t>{}));
#endif
            s /= s_norm;

            // if (node_id == 0) printf("s_norm = %f\n", s_norm);
            // printf("s = (%f, %f, %f)\n", s[0], s[1], s[2]);
        }
    }
};

template <ForcefieldType FFT, typename T, typename K>
SyclEvent compute_hessians(SyclQueue& Q, FullereneBatchView<T,K> B, Span<T> hess, Span<K> cols){
    TEMPLATE_TYPEDEFS(T,K);
    if (hess.size() < 90*B.N_*B.size() || cols.size() < 90*B.N_*B.size()) throw std::runtime_error("compute_hessians: hess and cols buffers must be of size >= 90*N*size");
    SyclEventImpl hessians_finished = Q->submit([&](sycl::handler& h){
        auto X_acc = B.d_.X_cubic_;
        auto cubic_neighbours_acc = B.d_.A_cubic_;

        sycl::local_accessor<coord3d, 1> X(B.N_,h);
        sycl::local_accessor<real_t, 1> sdata(3*B.N_,h);

        auto N = B.N_;
        auto size = B.size();
        h.parallel_for(sycl::nd_range(sycl::range{size*N}, sycl::range{N}), [=](sycl::nd_item<1> nditem){
            auto cta = nditem.get_group();
            auto tid = nditem.get_local_linear_id();
            auto bid = nditem.get_group_linear_id();
            if (!(B[bid].m_.flags_.get() & StatusFlag::FULLERENEGRAPH_PREPARED)) return; //Skip if cubic graph is not initialized

            Constants<T,K> constants (cubic_neighbours_acc, K(tid));
            NodeNeighbours nodeG(cubic_neighbours_acc, K(tid));
            X[tid] = X_acc[bid*N + tid];
            ForceField FF = ForceField<FFT,T,K>(nodeG, constants, cta, sdata.get_pointer());
            auto hessian = FF.hessian(Span<coord3d>(X.get_pointer(), N));
            int n_cols = 10*3;
            int n_rows = N*3;
            int hess_stride = n_cols*n_rows;
            int toff = bid*hess_stride + tid*n_cols*3;
            for (size_t i = 0; i < 3; i++) //rows
            for (size_t j = 0; j < 10; j++) //cols / 3
            {   
                for (size_t k = 0; k < 3; k++)
                {   
                    cols[toff + i*n_cols + j*3 + k] = hessian.indices[j]*3 + k;
                    hess[toff + i*n_cols + j*3 + k] = hessian.A[j][i][k];
                }    
            }
            sycl::group_barrier(cta);
            //Enforce symmetry

            for(int ii = tid; ii < n_cols*n_rows; ii += N){
                int i = ii / n_cols;
                int jj = ii % n_cols;
                int j = cols[bid*hess_stride + i*n_cols + jj];
                int ix = 0;
                
                if(i < j){
                    while (ix < n_cols && cols[bid*hess_stride + j*n_cols + ix] != i) {ix++;}
                    int jx = cols[bid*hess_stride + j*n_cols + ix];

                    T val = 0.5*(hess[bid*hess_stride + i*n_cols + jj] + hess[bid*hess_stride + j*n_cols + ix]);
                    hess[bid*hess_stride + i*n_cols + jj] = val;
                    hess[bid*hess_stride + j*n_cols + ix] = val;
                }
            }

        });
    });
    return SyclEvent(std::move(hessians_finished));
}

template <ForcefieldType FFT, typename T, typename K>
SyclEvent HessianFunctor<FFT,T,K>::compute(SyclQueue& Q, FullereneBatchView<T,K> B, Span<T> hess, Span<K> cols){
    return compute_hessians<FFT, T, K>(Q, B, hess, cols);
}

template <ForcefieldType FFT, typename T, typename K>
SyclEvent HessianFunctor<FFT,T,K>::compute(SyclQueue& Q, Fullerene<T,K> B, Span<T> hess, Span<K> cols, Span<K> indices){
    //TODO: Implement compute for single isomer
    throw std::logic_error("HessianFunctor::compute not implemented for single isomer");
}

template struct HessianFunctor<PEDERSEN, float, uint16_t>;
template struct HessianFunctor<PEDERSEN, double, uint16_t>;
template struct HessianFunctor<PEDERSEN, float, uint32_t>;
template struct HessianFunctor<PEDERSEN, double, uint32_t>;

/* template void compute_hessians<PEDERSEN, float, uint16_t>(sycl::queue& Q, IsomerBatch<float,uint16_t>& B, sycl::buffer<float,1>& hess, sycl::buffer<uint16_t,1>& cols, const LaunchPolicy policy);
template void compute_hessians<PEDERSEN, double, uint16_t>(sycl::queue& Q, IsomerBatch<double,uint16_t>& B, sycl::buffer<double,1>& hess, sycl::buffer<uint16_t,1>& cols, const LaunchPolicy policy); */
