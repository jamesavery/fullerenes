# Point Group 3D Representation Plan

## Goal

Add a method `Symmetry::representation_3d()` that returns the 3D rotation/reflection
matrices for each element of the automorphism group G. This enables genuine
symmetry-constrained optimization in the iDT embed_3d() pipeline.

## Data Structure

```cpp
struct Representation3D {
    vector<matrix3d> R;    // R[i] = 3D matrix for G[i]
    // Invariant: R[i] * R[j] == R[k] whenever G[i]*G[j] == G[k]
    // det(R[i]) = +1 for proper, -1 for improper rotations
};
```

## Algorithm

### Step 1: Generate standard 3D matrices from the point group type

For each of the 28 fullerene point groups, generate the full set of 3D matrices
from known generators using standard crystallographic conventions:

- Principal Cn axis along z
- For dihedral groups: one C2' axis along x
- For cubic/icosahedral: standard orientation (C3 along (1,1,1) for T/Td/Th,
  golden-ratio icosahedron for I/Ih)

Generator matrices:
- Rz(theta): rotation by theta about z
- C2x: diag(1,-1,-1) -- 180 deg about x
- sigma_h: diag(1,1,-1) -- reflection in xy-plane
- sigma_v: diag(1,-1,1) -- reflection in xz-plane
- inversion: diag(-1,-1,-1)
- sigma_d: reflection in the xz-plane rotated by pi/(2n) (bisects C2' axes)

Groups and their generators:

| Group   | |G| | Generators                                |
|---------|-----|-------------------------------------------|
| C1      |  1  | {E}                                       |
| C2      |  2  | Rz(pi)                                    |
| Ci      |  2  | inversion                                 |
| Cs      |  2  | sigma_h                                   |
| C3      |  3  | Rz(2pi/3)                                 |
| S4      |  4  | S4 = Rz(pi/2) * sigma_h                   |
| C2v     |  4  | Rz(pi), sigma_v                           |
| C2h     |  4  | Rz(pi), sigma_h                           |
| D2      |  4  | Rz(pi), C2x                               |
| S6      |  6  | S6 = Rz(pi/3) * sigma_h                   |
| C3v     |  6  | Rz(2pi/3), sigma_v                        |
| C3h     |  6  | Rz(2pi/3), sigma_h                        |
| D3      |  6  | Rz(2pi/3), C2x                            |
| D2h     |  8  | Rz(pi), C2x, sigma_h                      |
| D2d     |  8  | S4z, C2x                                  |
| D5      | 10  | Rz(2pi/5), C2x                            |
| D6      | 12  | Rz(pi/3), C2x                             |
| D3h     | 12  | Rz(2pi/3), C2x, sigma_h                   |
| D3d     | 12  | S6z, C2x                                  |
| T       | 12  | C3(1,1,1), C2z                            |
| D5h     | 20  | Rz(2pi/5), C2x, sigma_h                   |
| D5d     | 20  | S10z, C2x                                 |
| D6h     | 24  | Rz(pi/3), C2x, sigma_h                    |
| D6d     | 24  | S12z, C2x                                 |
| Td      | 24  | C3(1,1,1), S4z                            |
| Th      | 24  | C3(1,1,1), C2z, inversion                 |
| I       | 60  | C5(gold), C3(1,1,1)                       |
| Ih      | 120 | C5(gold), C3(1,1,1), inversion            |

### Step 2: Match permutations to matrices

Find an isomorphism phi: {3D matrices} -> {permutations} such that the
multiplication tables match.

Algorithm:
1. Classify each permutation by: (order, orientation_character).
   - order = pi.order()
   - orientation = reverses_orientation(pi) ? -1 : +1
2. Classify each matrix by: (order, det).
   - det = +1 proper, -1 improper
3. Group by (order, character) -- these are "type classes".
4. For each type class, try to match matrices to permutations.
   - Start with the identity (trivially matched).
   - For generators: pick a candidate permutation of matching type, assign it.
   - Generate the full matching by closure under multiplication.
   - Verify the multiplication table is consistent.
   - If inconsistent, backtrack and try another candidate.

For the 28 fullerene groups (max order 120), this is fast: the generator
matching has very few candidates per type class.

### Step 3: Validation

Verify the representation satisfies:
1. R[i] * R[j] == R[k] whenever G[i]*G[j] == G[k]  (homomorphism)
2. det(R[i]) == +1 iff G[i] preserves orientation   (character match)
3. |R[i]| == |G|  (same size)
4. For I/Ih: golden-ratio structure is a sanity check

## Files Modified

- include/fullerenes/symmetry.hh -- add Representation3D struct and method declaration
- src/c++/symmetry.cc -- implement representation_3d()
- tests/symmetry-test.cc -- add validation tests

## Testing

For each point group that occurs among fullerene isomers (at least C20-C60):
1. Compute Symmetry, call representation_3d()
2. Verify multiplication table consistency
3. Verify det matches orientation character
4. Verify all matrices are orthogonal (R^T R == I)
5. Verify group order matches

## Usage in embed_3d()

Once representation_3d() is available, the symmetry-constrained optimization
parameterizes only orbit representatives:

For each group element g with permutation pi_g and matrix R_g:
  x[pi_g(v)] = R_g * x[v]

This reduces 12*3 = 36 DOFs to (n_orbits * 3 - site_symmetry_constraints).
No MDS collapse is possible because orbit members are never placed independently.
