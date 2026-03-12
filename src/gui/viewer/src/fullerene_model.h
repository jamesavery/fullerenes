#pragma once

#include "mesh.h"
#include <vector>
#include <string>

class Deltahedron;

// Quality coloring modes
enum class ColorMode {
    Solid = 0,
    AngleError = 1,
    Convexity = 2,
    Degree = 3,
};

struct MeshData {
    std::vector<Vertex> vertices;
    std::vector<uint32_t> triangleIndices;
    std::vector<uint32_t> edgeIndices;
    glm::vec3 center;
    float radius;
};

// Symmetry validation result
struct SymmetryCheck {
    std::string pointGroup;       // Combinatorial point group (e.g. "Ih")
    int groupOrder = 0;           // |G| = number of automorphisms
    int geometricCount = 0;       // How many automorphisms are realized as isometries
    double maxRMSD = 0;           // Worst RMSD across all automorphisms
    double meanRMSD = 0;          // Mean RMSD across all automorphisms
    bool valid() const { return geometricCount == groupOrder; }
};

// Convert a Deltahedron into GPU-ready mesh data.
// quality_mode determines what per-vertex metric is stored in Vertex::quality.
MeshData deltahedronToMesh(const Deltahedron& D, ColorMode mode = ColorMode::Solid);

// Check if the 3D geometry respects the combinatorial symmetry group.
// For each automorphism g, checks that coords[g[i]] ≈ R·coords[i] + t for some isometry R,t.
SymmetryCheck checkSymmetry(const Deltahedron& D);
