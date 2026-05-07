#include "fullerenes/triangulation.hh"
#include "fullerenes/symmetry.hh"
#include "fullerenes/buckygen-wrapper.hh"
#include <sys/wait.h>
#include <csignal>

// Compare symmetry group order from buckygen graph vs spiral-reconstructed graph.
// The spiral-reconstructed graph has P0 = identity, so permutation_representation()
// is known to be correct for it (the original code worked for spiral-constructed graphs).
int test_symmetry_order(int N)
{
  signal(SIGTERM, SIG_IGN);

  Graph g;
  BuckyGen::buckygen_queue Q = BuckyGen::start(N, false, false);

  int count = 0, fails = 0;
  while(BuckyGen::next_fullerene(Q, g)){
    count++;
    Triangulation T(g);

    // Get canonical spiral from buckygen graph
    vector<int> spiral;
    jumplist_t jumps;
    if(!T.get_spiral(spiral, jumps)){ fails++; continue; }

    // Symmetry from buckygen graph (uses P0⁻¹∘P with validation)
    Symmetry S_bucky(T);
    int order_bucky = S_bucky.permutation_representation().size();

    // Symmetry from spiral-reconstructed graph (P0 = identity, known correct)
    Triangulation T2(spiral, jumps);
    Symmetry S_spiral(T2);
    int order_spiral = S_spiral.permutation_representation().size();

    if(order_bucky != order_spiral){
      fprintf(stderr, "  ORDER MISMATCH: N=%d #%d: buckygen=%d spiral=%d\n",
              N, count, order_bucky, order_spiral);
      fails++;
    }
  }

  fprintf(stderr, "N=%3d: %6d isomers, %d order mismatches\n", N, count, fails);
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
      int fails = test_symmetry_order(N);
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
  fprintf(stderr, "\nTotal: %d order mismatches\n", total_failures);
  return total_failures > 0 ? 1 : 0;
}
