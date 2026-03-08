#include "fullerenes/triangulation.hh"
#include "fullerenes/buckygen-wrapper.hh"
#include <sys/wait.h>
#include <csignal>

int test_triple(const Triangulation& T, node_t u, node_t v, node_t w, bool general)
{
  vector<int> spiral;
  jumplist_t jumps;
  vector<node_t> perm;
  bool ok = T.get_spiral_implementation(u, v, w, spiral, jumps, perm, general);
  return ok ? 0 : 1;
}

void diagnose_isomer(const Triangulation& T, int N, int idx)
{
  fprintf(stderr, "\n  Diagnosing N=%d isomer #%d (N_tri=%d, N_vert=%zu):\n",
          N, idx, T.N, T.triangles().size());

  // Print degree sequence
  fprintf(stderr, "  Degrees:");
  for(node_t u = 0; u < T.N; u++)
    if(T.neighbours[u].size() != 6) fprintf(stderr, " %d(%zu)", u, T.neighbours[u].size());
  fprintf(stderr, "\n");

  // Try every starting triple individually in a forked process
  int triple_count = 0, crash_count = 0;
  for(node_t u = 0; u < T.N; u++){
    for(node_t v : T.neighbours[u]){
      node_t w_prev = T.prev(v, u);
      node_t w_next = T.next(v, u);

      for(int k = 0; k < 2; k++){
        node_t w = (k == 0) ? w_prev : w_next;
        triple_count++;

        for(int gen = 0; gen < 2; gen++){
          pid_t pid = fork();
          if(pid == 0){
            int r = test_triple(T, u, v, w, gen);
            _exit(r);
          }
          int status;
          waitpid(pid, &status, 0);
          if(WIFSIGNALED(status)){
            fprintf(stderr, "    CRASH: triple (%d,%d,%d) general=%s signal=%d\n",
                    u, v, w, gen ? "true" : "false", WTERMSIG(status));
            crash_count++;
          }
        }
      }
    }
  }
  fprintf(stderr, "  Tested %d triples (x2 general), %d crashes\n", triple_count, crash_count);
}

int main(int argc, char** argv)
{
  int N = argc > 1 ? atoi(argv[1]) : 44;
  int target_isomer = argc > 2 ? atoi(argv[2]) : 0; // 0 = test all
  fprintf(stderr, "Diagnosing canonical spiral crashes for N=%d\n", N);

  signal(SIGTERM, SIG_IGN);

  Graph g;
  BuckyGen::buckygen_queue Q = BuckyGen::start(N, false, false);

  int count = 0;
  while(BuckyGen::next_fullerene(Q, g)){
    count++;
    if(target_isomer > 0 && count != target_isomer) continue;

    Triangulation T(g);

    // Try canonical search in a forked process
    pid_t pid = fork();
    if(pid == 0){
      vector<int> spiral;
      jumplist_t jumps;
      bool ok = T.get_spiral(spiral, jumps);
      _exit(ok ? 0 : 1);
    }
    int status;
    waitpid(pid, &status, 0);
    if(WIFSIGNALED(status)){
      fprintf(stderr, "N=%d isomer #%d: CANONICAL CRASHED (signal %d)\n",
              N, count, WTERMSIG(status));
      diagnose_isomer(T, N, count);
    }

    if(target_isomer > 0) break;
  }
  fprintf(stderr, "\nDone: tested %d isomers for N=%d\n", count, N);
  return 0;
}
