#include "fullerene_model.h"
#include <fullerenes/deltahedron.hh>
#include <fullerenes/symmetry.hh>
#include <fullerenes/geometry.hh>
#include <glm/glm.hpp>
#include <cmath>
#include <algorithm>
#include <numeric>

MeshData deltahedronToMesh(const Deltahedron& D, ColorMode mode) {
    MeshData data;
    int N = D.N;

    // Convert positions (double -> float) and compute bounding sphere
    glm::vec3 center(0);
    for (int i = 0; i < N; i++) {
        center += glm::vec3((float)D.points[i][0], (float)D.points[i][1], (float)D.points[i][2]);
    }
    center /= (float)N;
    data.center = center;

    float maxR = 0;
    for (int i = 0; i < N; i++) {
        glm::vec3 p((float)D.points[i][0], (float)D.points[i][1], (float)D.points[i][2]);
        maxR = std::max(maxR, glm::length(p - center));
    }
    data.radius = maxR;

    // Compute per-vertex normals (area-weighted face normals)
    std::vector<glm::vec3> normals(N, glm::vec3(0));
    for (const auto& tri : D.triangles) {
        glm::vec3 v0((float)D.points[tri[0]][0], (float)D.points[tri[0]][1], (float)D.points[tri[0]][2]);
        glm::vec3 v1((float)D.points[tri[1]][0], (float)D.points[tri[1]][1], (float)D.points[tri[1]][2]);
        glm::vec3 v2((float)D.points[tri[2]][0], (float)D.points[tri[2]][1], (float)D.points[tri[2]][2]);
        glm::vec3 faceNormal = glm::cross(v1 - v0, v2 - v0);  // area-weighted (not normalized)
        normals[tri[0]] += faceNormal;
        normals[tri[1]] += faceNormal;
        normals[tri[2]] += faceNormal;
    }
    for (int i = 0; i < N; i++) {
        float len = glm::length(normals[i]);
        if (len > 1e-10f) normals[i] /= len;
    }

    // Compute per-vertex quality metric
    std::vector<float> quality(N, 0.0f);
    if (mode == ColorMode::AngleError) {
        // Max angle error at this vertex across all incident triangles
        for (const auto& tri : D.triangles) {
            for (int c = 0; c < 3; c++) {
                int vi = tri[c];
                coord3d va = D.points[tri[(c+1)%3]] - D.points[tri[c]];
                coord3d vb = D.points[tri[(c+2)%3]] - D.points[tri[c]];
                double angle = coord3d::angle(va, vb) * 180.0 / M_PI;
                float err = (float)(std::abs(angle - 60.0) / 60.0);
                quality[vi] = std::max(quality[vi], err);
            }
        }
        // Absolute scale: 0% error = 0 (blue), 1%+ error = 1 (red)
        for (auto& q : quality) q = std::clamp(q / 0.01f, 0.0f, 1.0f);
    } else if (mode == ColorMode::Convexity) {
        // Signed height h (negative = concave)
        for (int v = 0; v < N; v++) {
            int d = (int)D.neighbours[v].size();
            coord3d cen(0,0,0);
            for (int i = 0; i < d; i++) cen += D.points[D.neighbours[v][i]];
            cen /= (double)d;
            coord3d nf(0,0,0);
            for (int i = 0; i < d; i++) {
                coord3d e1 = D.points[D.neighbours[v][i]] - D.points[v];
                coord3d e2 = D.points[D.neighbours[v][(i+1)%d]] - D.points[v];
                nf += e1.cross(e2);
            }
            double nl = nf.norm();
            if (nl > 1e-15) {
                double h = (D.points[v] - cen).dot(nf / nl);
                // Map: h < 0 (concave) -> 1.0 (red), h > 0 (convex) -> 0.0 (blue)
                quality[v] = std::clamp((float)(-h * 50.0), 0.0f, 1.0f);
            }
        }
    } else if (mode == ColorMode::Degree) {
        for (int v = 0; v < N; v++) {
            quality[v] = ((int)D.neighbours[v].size() == 5) ? 1.0f : 0.0f;
        }
    }

    // Build vertex buffer
    data.vertices.resize(N);
    for (int i = 0; i < N; i++) {
        data.vertices[i].pos = glm::vec3((float)D.points[i][0], (float)D.points[i][1], (float)D.points[i][2]);
        data.vertices[i].normal = normals[i];
        data.vertices[i].quality = quality[i];
    }

    // Triangle indices
    data.triangleIndices.reserve(D.triangles.size() * 3);
    for (const auto& tri : D.triangles) {
        data.triangleIndices.push_back((uint32_t)tri[0]);
        data.triangleIndices.push_back((uint32_t)tri[1]);
        data.triangleIndices.push_back((uint32_t)tri[2]);
    }

    // Edge indices (line list, each edge once)
    for (int u = 0; u < N; u++) {
        for (int v : D.neighbours[u]) {
            if (v > u) {
                data.edgeIndices.push_back((uint32_t)u);
                data.edgeIndices.push_back((uint32_t)v);
            }
        }
    }

    return data;
}

// LAPACK SVD for 3x3 matrix
extern "C" void dgesvd_(const char* jobu, const char* jobvt, const int* m, const int* n,
                         double* A, const int* lda, double* S, double* U, const int* ldu,
                         double* VT, const int* ldvt, double* work, const int* lwork, int* info);

// Compute RMSD between point set P and permuted set Q[i] = P[perm[i]]
// after optimal rigid alignment (Kabsch algorithm).
static double kabschRMSD(const std::vector<coord3d>& P, const std::vector<int>& perm) {
    int N = (int)P.size();
    if (N == 0) return 0;

    // Compute centroids
    coord3d cP(0,0,0), cQ(0,0,0);
    for (int i = 0; i < N; i++) {
        cP += P[i];
        cQ += P[perm[i]];
    }
    cP = cP / (double)N;
    cQ = cQ / (double)N;

    // Build 3x3 cross-covariance matrix H = sum (P[i]-cP) * (Q[i]-cQ)^T
    // H[r][c] = sum_i (P[i][r] - cP[r]) * (Q[i][c] - cQ[c])
    double H[9] = {};  // column-major for LAPACK
    for (int i = 0; i < N; i++) {
        coord3d p = P[i] - cP;
        coord3d q = P[perm[i]] - cQ;
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 3; c++)
                H[c*3 + r] += p[r] * q[c];  // column-major
    }

    // SVD: H = U * S * VT
    int m = 3, n = 3, lda = 3, ldu = 3, ldvt = 3, info;
    double S[3], U[9], VT[9];
    double work_query;
    int lwork = -1;
    dgesvd_("A", "A", &m, &n, H, &lda, S, U, &ldu, VT, &ldvt, &work_query, &lwork, &info);
    lwork = (int)work_query;
    std::vector<double> work(lwork);
    dgesvd_("A", "A", &m, &n, H, &lda, S, U, &ldu, VT, &ldvt, work.data(), &lwork, &info);

    if (info != 0) return 1e10;

    // R = V * U^T, ensuring det(R) > 0 (proper rotation)
    // det(V * U^T) = det(V) * det(U^T) = det(V) * det(U)
    // If negative, flip sign of last column of U
    // det of 3x3 column-major matrix M: M[0]*M[4]*M[8] + ...
    auto det3 = [](const double* M) {
        return M[0]*(M[4]*M[8]-M[7]*M[5]) - M[3]*(M[1]*M[8]-M[7]*M[2]) + M[6]*(M[1]*M[5]-M[4]*M[2]);
    };
    double d = det3(U) * det3(VT);
    if (d < 0) {
        // Flip last column of U
        U[6] = -U[6]; U[7] = -U[7]; U[8] = -U[8];
    }

    // R = V * U^T.  V = VT^T, so R = VT^T * U^T = (U * VT)^T
    // R[r][c] = sum_k V[r][k] * U[c][k] = sum_k VT[k*3+r] * U[k*3+c]
    double R[9];  // column-major
    for (int r = 0; r < 3; r++)
        for (int c = 0; c < 3; c++) {
            double sum = 0;
            for (int k = 0; k < 3; k++)
                sum += VT[k*3 + r] * U[k*3 + c];
            R[c*3 + r] = sum;
        }

    // Compute RMSD: sqrt( sum_i ||R*(P[i]-cP) - (Q[i]-cQ)||^2 / N )
    double sse = 0;
    for (int i = 0; i < N; i++) {
        coord3d p = P[i] - cP;
        coord3d q = P[perm[i]] - cQ;
        // Rp = R * p (column-major)
        coord3d Rp;
        for (int r = 0; r < 3; r++) {
            Rp[r] = 0;
            for (int c = 0; c < 3; c++)
                Rp[r] += R[c*3 + r] * p[c];
        }
        coord3d diff = Rp - q;
        sse += diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2];
    }
    return std::sqrt(sse / N);
}

SymmetryCheck checkSymmetry(const Deltahedron& D) {
    SymmetryCheck result;

    // Compute combinatorial symmetry from the triangulation
    Symmetry sym(static_cast<const Triangulation&>(D));

    result.pointGroup = sym.point_group().to_string();
    result.groupOrder = (int)sym.G.size();

    // For each automorphism, check if it's realized as a geometric isometry
    double maxRMSD = 0, sumRMSD = 0;
    int nGeometric = 0;
    double tol = 1e-6;  // RMSD tolerance for "isometry"

    // Scale tolerance by mean edge length
    double meanL = 0;
    int ne = 0;
    for (int u = 0; u < D.N; u++)
        for (int v : D.neighbours[u])
            if (v > u) { meanL += (D.points[u] - D.points[v]).norm(); ne++; }
    if (ne > 0) meanL /= ne;
    tol = meanL * 1e-6;

    for (const auto& g : sym.G) {
        double rmsd = kabschRMSD(D.points, g);
        maxRMSD = std::max(maxRMSD, rmsd);
        sumRMSD += rmsd;
        if (rmsd < tol) nGeometric++;
    }

    result.geometricCount = nGeometric;
    result.maxRMSD = maxRMSD;
    result.meanRMSD = result.groupOrder > 0 ? sumRMSD / result.groupOrder : 0;

    return result;
}
