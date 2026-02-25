#include "fullerenes/triangulation.hh"
#include "fullerenes/symmetry.hh"
#include "fullerenes/buckygen-wrapper.hh"
#include <sys/wait.h>
#include <csignal>

// Test 1: Single-triple spiral extraction
int test_single_spiral(int N)
{
  signal(SIGTERM, SIG_IGN);

  Graph g;
  BuckyGen::buckygen_queue Q = BuckyGen::start(N, false, false);

  int count = 0, fails = 0;
  while(BuckyGen::next_fullerene(Q, g)){
    count++;
    Triangulation T(g);

    int expected_tris = 2 * (T.N - 2);
    if((int)T.triangles.size() != expected_tris){
      fprintf(stderr, "  TRI FAIL: N=%d #%d: got %zu expected %d\n",
              N, count, T.triangles.size(), expected_tris);
      fails++;
      continue;
    }

    node_t u = 0;
    node_t v = T.neighbours[u][0];
    node_t w = T.next(u, v);
    vector<int> spiral;
    jumplist_t jumps;
    vector<node_t> perm;
    bool ok = T.get_spiral(u, v, w, spiral, jumps, perm);
    if(!ok){
      fprintf(stderr, "  SINGLE SPIRAL FAIL: N=%d #%d\n", N, count);
      fails++;
    }
  }

  fprintf(stderr, "  single-triple: N=%3d: %6d isomers, %d fails\n", N, count, fails);
  return fails;
}

// Test 2: Canonical spiral search (tries ALL starting triples)
int test_canonical_spiral(int N)
{
  signal(SIGTERM, SIG_IGN);

  Graph g;
  BuckyGen::buckygen_queue Q = BuckyGen::start(N, false, false);

  int count = 0, fails = 0;
  while(BuckyGen::next_fullerene(Q, g)){
    count++;
    Triangulation T(g);

    vector<int> spiral;
    jumplist_t jumps;
    bool ok = T.get_spiral(spiral, jumps);
    if(!ok){
      fprintf(stderr, "  CANONICAL FAIL: N=%d #%d\n", N, count);
      fails++;
    }
  }

  fprintf(stderr, "  canonical:     N=%3d: %6d isomers, %d fails\n", N, count, fails);
  return fails;
}

// Test 3: Spiral round-trip (extract spiral, wind up into Triangulation, extract again, compare)
int test_roundtrip(int N)
{
  signal(SIGTERM, SIG_IGN);

  Graph g;
  BuckyGen::buckygen_queue Q = BuckyGen::start(N, false, false);

  int count = 0, fails = 0;
  while(BuckyGen::next_fullerene(Q, g)){
    count++;
    Triangulation T(g);

    // Extract canonical spiral
    vector<int> spiral1;
    jumplist_t jumps1;
    bool ok1 = T.get_spiral(spiral1, jumps1);
    if(!ok1){
      fails++;
      continue;
    }

    // Wind up into new Triangulation
    Triangulation T2(spiral1, jumps1);

    // Extract canonical spiral from the reconstructed triangulation
    vector<int> spiral2;
    jumplist_t jumps2;
    bool ok2 = T2.get_spiral(spiral2, jumps2);
    if(!ok2){
      fprintf(stderr, "  ROUNDTRIP FAIL (2nd extract): N=%d #%d\n", N, count);
      fails++;
      continue;
    }

    // Compare
    if(spiral1 != spiral2 || jumps1 != jumps2){
      fprintf(stderr, "  ROUNDTRIP MISMATCH: N=%d #%d\n", N, count);
      fails++;
    }
  }

  fprintf(stderr, "  roundtrip:     N=%3d: %6d isomers, %d fails\n", N, count, fails);
  return fails;
}

// Test 4: Symmetry detection (point group from buckygen graph)
int test_symmetry(int N)
{
  signal(SIGTERM, SIG_IGN);

  Graph g;
  BuckyGen::buckygen_queue Q = BuckyGen::start(N, false, false);

  int count = 0, fails = 0;
  while(BuckyGen::next_fullerene(Q, g)){
    count++;
    Triangulation T(g);

    Symmetry S(T);
    auto pg = S.point_group();

    // Validate: point group must have a recognized type
    if(pg.sym_type == PointGroup::UNKNOWN){
      fprintf(stderr, "  SYMMETRY FAIL (unknown group): N=%d #%d\n", N, count);
      fails++;
      continue;
    }

    // Verify all automorphisms are valid permutations
    auto perms = S.permutation_representation();
    for(size_t i = 0; i < perms.size(); i++){
      // Check permutation is valid (maps to valid nodes)
      bool valid = true;
      vector<bool> seen(T.N, false);
      for(int j = 0; j < T.N; j++){
        int p = perms[i][j];
        if(p < 0 || p >= T.N || seen[p]){
          valid = false;
          break;
        }
        seen[p] = true;
      }
      if(!valid){
        fprintf(stderr, "  SYMMETRY FAIL (invalid perm %zu): N=%d #%d\n", i, N, count);
        fails++;
        break;
      }

      // Check permutation preserves adjacency
      for(int u = 0; u < T.N; u++){
        for(node_t v : T.neighbours[u]){
          int pu = perms[i][u], pv = perms[i][v];
          if(T.arc_ix(pu, pv) < 0){
            fprintf(stderr, "  SYMMETRY FAIL (perm %zu breaks adj %d->%d => %d->%d): N=%d #%d\n",
                    i, u, v, pu, pv, N, count);
            fails++;
            goto next_isomer;
          }
        }
      }
    }
    next_isomer:;
  }

  fprintf(stderr, "  symmetry:      N=%3d: %6d isomers, %d fails\n", N, count, fails);
  return fails;
}

int run_forked(const char* test_name, int N, int (*test_fn)(int))
{
  pid_t pid = fork();
  if(pid == 0){
    int fails = test_fn(N);
    _exit(fails > 127 ? 127 : fails);
  }
  int status;
  waitpid(pid, &status, 0);
  if(WIFSIGNALED(status)){
    fprintf(stderr, "  %-14s N=%3d: CRASHED (signal %d)\n", test_name, N, WTERMSIG(status));
    return 1;
  } else if(WIFEXITED(status)){
    return WEXITSTATUS(status);
  }
  return 1;
}

int main()
{
  int sizes[] = {20, 24, 26, 28, 30, 32, 34, 36, 38, 40,
                 42, 44, 46, 48, 50, 52, 54, 56, 58, 60};

  fprintf(stderr, "=== Test 1: Single-triple spiral extraction ===\n");
  int single_fails = 0;
  for(int N : sizes)
    single_fails += run_forked("single-triple", N, test_single_spiral);
  fprintf(stderr, "Single-triple total: %d failures\n\n", single_fails);

  fprintf(stderr, "=== Test 2: Canonical spiral search ===\n");
  int canonical_fails = 0;
  for(int N : sizes)
    canonical_fails += run_forked("canonical", N, test_canonical_spiral);
  fprintf(stderr, "Canonical total: %d failures\n\n", canonical_fails);

  fprintf(stderr, "=== Test 3: Spiral round-trip ===\n");
  int roundtrip_fails = 0;
  for(int N : sizes)
    roundtrip_fails += run_forked("roundtrip", N, test_roundtrip);
  fprintf(stderr, "Round-trip total: %d failures\n\n", roundtrip_fails);

  fprintf(stderr, "=== Test 4: Symmetry detection ===\n");
  int symmetry_fails = 0;
  for(int N : sizes)
    symmetry_fails += run_forked("symmetry", N, test_symmetry);
  fprintf(stderr, "Symmetry total: %d failures\n\n", symmetry_fails);

  fprintf(stderr, "=== Summary ===\n");
  fprintf(stderr, "Single-triple: %d failures\n", single_fails);
  fprintf(stderr, "Canonical:     %d failures\n", canonical_fails);
  fprintf(stderr, "Round-trip:    %d failures\n", roundtrip_fails);
  fprintf(stderr, "Symmetry:      %d failures\n", symmetry_fails);

  int total = single_fails + canonical_fails + roundtrip_fails + symmetry_fails;
  return total > 0 ? 1 : 0;
}
