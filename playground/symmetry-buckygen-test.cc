#include "fullerenes/triangulation.hh"
#include "fullerenes/symmetry.hh"
#include "fullerenes/buckygen-wrapper.hh"
#include <sys/wait.h>
#include <unistd.h>
#include <csignal>

// All tests in a single pass per isomer.
int test_all(int N)
{
  signal(SIGTERM, SIG_IGN);

  Graph g;
  BuckyGen::buckygen_queue Q = BuckyGen::start(N, false, false);

  int count = 0, fails = 0;
  while(BuckyGen::next_fullerene(Q, g)){
    count++;
    Triangulation T(g);

    // --- Triangulation sanity ---
    int expected_tris = 2 * (T.N - 2);
    if((int)T.triangles().size() != expected_tris){
      fprintf(stderr, "  TRI FAIL: N=%d #%d: got %zu expected %d\n",
              N, count, T.triangles().size(), expected_tris);
      fails++;
      continue;
    }

    // --- Single-triple spiral extraction ---
    {
      node_t u = 0, v = T[u][0], w = T.next(u, v);
      vector<int> spiral;
      jumplist_t jumps;
      vector<node_t> perm;
      if(!T.get_spiral(u, v, w, spiral, jumps, perm)){
        fprintf(stderr, "  SINGLE SPIRAL FAIL: N=%d #%d\n", N, count);
        fails++;
        continue;
      }
    }

    // --- Canonical spiral search ---
    vector<int> spiral1;
    jumplist_t jumps1;
    if(!T.get_spiral(spiral1, jumps1)){
      fprintf(stderr, "  CANONICAL FAIL: N=%d #%d\n", N, count);
      fails++;
      continue;
    }

    // --- Spiral round-trip ---
    {
      Triangulation T2(spiral1, jumps1);
      vector<int> spiral2;
      jumplist_t jumps2;
      if(!T2.get_spiral(spiral2, jumps2)){
        fprintf(stderr, "  ROUNDTRIP FAIL (2nd extract): N=%d #%d\n", N, count);
        fails++;
        continue;
      }
      if(spiral1 != spiral2 || jumps1 != jumps2){
        fprintf(stderr, "  ROUNDTRIP MISMATCH: N=%d #%d\n", N, count);
        fails++;
        continue;
      }
    }

    // --- Symmetry detection ---
    {
      Symmetry S(T);
      auto pg = S.point_group();
      if(pg.sym_type == PointGroup::UNKNOWN){
        fprintf(stderr, "  SYMMETRY FAIL (unknown group): N=%d #%d\n", N, count);
        fails++;
        continue;
      }

      auto perms = S.permutation_representation();
      for(size_t i = 0; i < perms.size(); i++){
        // Check permutation is a valid bijection
        vector<bool> seen(T.N, false);
        bool valid = true;
        for(int j = 0; j < T.N; j++){
          int p = perms[i][j];
          if(p < 0 || p >= T.N || seen[p]){ valid = false; break; }
          seen[p] = true;
        }
        if(!valid){
          fprintf(stderr, "  SYMMETRY FAIL (invalid perm %zu): N=%d #%d\n", i, N, count);
          fails++;
          goto next_isomer;
        }

        // Check permutation preserves adjacency
        for(int u = 0; u < T.N; u++)
          for(node_t v : T[u])
            if(T.arc_ix(perms[i][u], perms[i][v]) < 0){
              fprintf(stderr, "  SYMMETRY FAIL (perm %zu breaks adj %d->%d => %d->%d): N=%d #%d\n",
                      i, u, v, perms[i][u], perms[i][v], N, count);
              fails++;
              goto next_isomer;
            }
      }
    }
    next_isomer:;
  }

  fprintf(stderr, "N=%3d: %6d isomers, %d fails\n", N, count, fails);
  return fails;
}

int main()
{
  int sizes[] = {20, 24, 26, 28, 30, 32, 34, 36, 38, 40,
                 42, 44, 46, 48, 50, 52, 54, 56, 58, 60,
                 62, 64, 66, 68, 70, 72, 74, 76, 78, 80};

  int total_failures = 0;
  for(int N : sizes){
    pid_t pid = fork();
    if(pid == 0){
      int fails = test_all(N);
      _exit(fails > 127 ? 127 : fails);
    }
    int status;
    waitpid(pid, &status, 0);
    if(WIFSIGNALED(status)){
      fprintf(stderr, "N=%3d: CRASHED (signal %d)\n", N, WTERMSIG(status));
      total_failures++;
    } else if(WIFEXITED(status)){
      total_failures += WEXITSTATUS(status);
    }
  }
  fprintf(stderr, "\nTotal: %d failures\n", total_failures);
  return total_failures > 0 ? 1 : 0;
}
