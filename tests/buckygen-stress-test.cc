// buckygen-stress-test -- stress the BuckyGen child-process teardown.
//
// Exercises every failure mode the teardown rework addresses, adversarially and
// under load, and measures process/IPC hygiene. Each scenario asserts (CHECK ->
// exit 1 on failure); a hard SIGALRM watchdog (exit 2) guarantees no scenario
// can hang CI. Counts/sizes are CLI flags (never env vars) so the same binary is
// both a fast CTest smoke (defaults) and a heavy manual stress.
//
//   buckygen-stress-test [--size N] [--cycles K] [--reads R] [--queues M]
//                        [--herd-chunks C] [--herd-workers W] [--budget SECONDS]
//
// Scenarios:
//   churn            K start->read R->stop cycles; no zombie/queue/fd leak after.
//   concurrent       M live queues; close half; survivors must keep producing.
//   sigterm-handler  with a host SIGTERM handler installed, stop() must not hang.
//   double-stop      drain fully (internal stop) then stop() again: safe no-op.
//   herd             Nchunks >> Nworkers driven to completion: correct count and
//                    no zombie accumulation across worker recycling.
//   parent-death     a sub-parent dies without stop(); its buckygen grandchild
//                    must self-exit via the getppid poll (no orphan).

#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/msg.h>
#include <sys/wait.h>
#include <unistd.h>
#include <signal.h>
#include <time.h>
#include <dirent.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "fullerenes/graph.hh"
#include "fullerenes/buckygen-wrapper.hh"

using namespace BuckyGen;

static const char* g_scenario = "init";   // named in the watchdog message

static void on_timeout(int) {
  const char* m = "\nbuckygen-stress: WATCHDOG TIMEOUT in scenario: ";
  write(2, m, strlen(m));
  write(2, g_scenario, strlen(g_scenario));
  write(2, " (a stop()/close() hung)\n", 25);
  _exit(2);
}

#define CHECK(cond, ...) do { if(!(cond)) {                                   \
    fprintf(stderr, "buckygen-stress: FAIL [%s]: ", g_scenario);             \
    fprintf(stderr, __VA_ARGS__); fprintf(stderr, "\n"); _exit(1);           \
  } } while(0)

// Reap any of our children that have exited; returns how many were collected.
// After a clean teardown this must be 0 (stop()/stop_all() already reaped them).
static int sweep_reap() {
  int n = 0;
  while (waitpid(-1, nullptr, WNOHANG) > 0) n++;
  return n;
}

// Count our children currently in the zombie ('Z') state, WITHOUT reaping them
// (so it observes accumulation rather than masking it). Linux /proc; 0 elsewhere.
static int count_zombie_children() {
#if defined(__linux__)
  DIR* d = opendir("/proc");
  if (!d) return -1;
  const int mypid = getpid();
  int n = 0;
  for (struct dirent* e; (e = readdir(d)); ) {
    if (e->d_name[0] < '0' || e->d_name[0] > '9') continue;
    char path[300];
    snprintf(path, sizeof path, "/proc/%s/stat", e->d_name);
    FILE* f = fopen(path, "r");
    if (!f) continue;
    char buf[512];
    if (fgets(buf, sizeof buf, f)) {
      // "pid (comm) state ppid ..." -- comm may contain spaces/parens, so parse
      // after the LAST ')'.
      char* rp = strrchr(buf, ')');
      char st; int ppid;
      if (rp && sscanf(rp + 1, " %c %d", &st, &ppid) == 2 && st == 'Z' && ppid == mypid)
        n++;
    }
    fclose(f);
  }
  closedir(d);
  return n;
#else
  return 0;
#endif
}

// Number of open file descriptors (Linux), for the leak check; -1 elsewhere.
static int count_open_fds() {
#if defined(__linux__)
  DIR* d = opendir("/proc/self/fd");
  if (!d) return -1;
  int n = 0;
  for (struct dirent* e; (e = readdir(d)); )
    if (e->d_name[0] != '.') n++;
  closedir(d);
  return n;
#else
  return -1;
#endif
}

// ---- scenarios -----------------------------------------------------------

static void scenario_churn(int N, int cycles, int reads) {
  g_scenario = "churn";
  const int fd0 = count_open_fds();
  for (int c = 0; c < cycles; c++) {
    buckygen_queue Q = start(N, false);
    Graph G;
    for (int i = 0; i < reads && next_fullerene(Q, G); i++) { /* read a few */ }
    stop(Q);
    CHECK(Q.pid == -1 && Q.qid == -1, "stop() left Q live (pid=%d qid=%d): not idempotent",
          (int)Q.pid, Q.qid);
  }
  const int leaked = sweep_reap();
  CHECK(leaked == 0, "%d unreaped child(ren) after churn (zombie leak)", leaked);
  const int fd1 = count_open_fds();
  CHECK(fd0 < 0 || fd1 <= fd0 + 2, "fd leak: %d -> %d across %d cycles", fd0, fd1, cycles);
  printf("  churn:           %d cycles x %d reads -- 0 leaked children, fds %d->%d\n",
         cycles, reads, fd0, fd1);
}

static void scenario_concurrent(int N, int M) {
  g_scenario = "concurrent";
  std::vector<buckygen_queue> qs;
  for (int i = 0; i < M; i++) qs.push_back(start(N, false));
  Graph G;
  for (int i = 0; i < M; i++)
    CHECK(next_fullerene(qs[i], G), "queue %d produced nothing", i);
  for (int i = 0; i < M; i += 2) stop(qs[i]);                 // close the even queues
  for (int i = 1; i < M; i += 2)                              // odd ones must survive
    CHECK(next_fullerene(qs[i], G), "survivor queue %d died when a sibling was closed", i);
  for (int i = 1; i < M; i += 2) stop(qs[i]);
  const int leaked = sweep_reap();
  CHECK(leaked == 0, "%d unreaped after concurrent", leaked);
  printf("  concurrent:      %d queues, closed %d, survivors kept producing\n", M, (M + 1) / 2);
}

static void scenario_sigterm_handler(int N, int cycles) {
  g_scenario = "sigterm-handler";
  struct sigaction sa{}, old{};
  sa.sa_handler = [](int) { /* a no-op host SIGTERM handler */ };
  sigaction(SIGTERM, &sa, &old);
  for (int c = 0; c < cycles; c++) {
    buckygen_queue Q = start(N, false);
    Graph G;
    next_fullerene(Q, G);
    stop(Q);                            // hangs (-> watchdog) if the child ran the inherited handler
  }
  sigaction(SIGTERM, &old, nullptr);
  const int leaked = sweep_reap();
  CHECK(leaked == 0, "%d unreaped after sigterm-handler", leaked);
  printf("  sigterm-handler: %d cycles, stop() never hung with a host handler installed\n", cycles);
}

static void scenario_double_stop(int N) {
  g_scenario = "double-stop";
  buckygen_queue Q = start(N, false);   // small N -> drains fast
  Graph G;
  int got = 0;
  while (next_fullerene(Q, G)) got++;   // last WORKER_FINISHED -> internal stop() nulls Q
  CHECK(Q.pid == -1, "Q not nulled after full drain (idempotency precondition)");
  stop(Q);                              // explicit second stop -> must be a safe no-op
  stop(Q);                              // and a third
  const int leaked = sweep_reap();
  CHECK(leaked == 0, "%d unreaped after double-stop", leaked);
  printf("  double-stop:     drained %d graph(s); repeated stop() is a safe no-op\n", got);
}

static int count_single(int N) {
  buckygen_queue Q = start(N, false);
  Graph G;
  int n = 0;
  while (next_fullerene(Q, G)) n++;
  return n;
}

static void scenario_herd(int N, int nchunks, int nworkers) {
  g_scenario = "herd";
  const int expected = count_single(N);
  sweep_reap();
  int peak_zombies = 0, got = 0;
  {
    buckyherd_queue HQ((size_t)N, (size_t)nchunks, (size_t)nworkers, false, false);
    Graph G;
    while (HQ.next_fullerene(G)) {
      got++;
      if ((got & 63) == 0) {                       // sample periodically
        int z = count_zombie_children();
        if (z > peak_zombies) peak_zombies = z;
      }
    }
  }                                                 // HQ dtor -> stop_all()
  const int leaked = sweep_reap();
  CHECK(got == expected, "herd count %d != single-queue count %d", got, expected);
  CHECK(peak_zombies <= nworkers + 1,
        "zombie accumulation across recycling: peak=%d (> %d); workers not reaped on recycle",
        peak_zombies, nworkers + 1);
  CHECK(leaked == 0, "%d unreaped after herd", leaked);
  printf("  herd:            %d chunks x %d workers -> %d graphs (== single-queue), peak zombies=%d\n",
         nchunks, nworkers, got, peak_zombies);
}

static void scenario_parent_death(int N, int budget_s) {
  g_scenario = "parent-death";
  struct { pid_t pid; int qid; } info{};
  int pfd[2];
  CHECK(pipe(pfd) == 0, "pipe failed");
  pid_t sub = fork();
  CHECK(sub >= 0, "fork failed");
  if (sub == 0) {                                   // sub-parent
    close(pfd[0]);
    buckygen_queue Q = start(N, false);
    Graph G;
    next_fullerene(Q, G);                           // read ONE graph, then leave the rest
    info.pid = Q.pid; info.qid = Q.qid;
    ssize_t w = write(pfd[1], &info, sizeof info);
    (void)w;
    _exit(0);                                        // DIE without stop() -> grandchild orphaned
  }
  close(pfd[1]);
  CHECK(read(pfd[0], &info, sizeof info) == (ssize_t)sizeof info, "no grandchild info from sub-parent");
  close(pfd[0]);
  waitpid(sub, nullptr, 0);                          // reap the sub-parent; grandchild now orphaned

  bool gone = false;                                 // grandchild must self-exit (getppid poll)
  for (int i = 0; i < budget_s * 20; i++) {
    if (kill(info.pid, 0) != 0) { gone = true; break; }   // ESRCH -> it exited
    struct timespec ts{0, 50 * 1000 * 1000};         // 50 ms
    nanosleep(&ts, nullptr);
  }
  msgctl(info.qid, IPC_RMID, 0);                     // clean the queue the dead sub-parent leaked
  CHECK(gone, "orphaned buckygen child %d did NOT self-exit; getppid poll failed", (int)info.pid);
  printf("  parent-death:    orphaned grandchild self-exited via the getppid poll\n");
}

// ---- driver --------------------------------------------------------------

int main(int argc, char** argv) {
  setvbuf(stdout, nullptr, _IONBF, 0);               // unbuffered, incremental progress
  setvbuf(stderr, nullptr, _IONBF, 0);

  int N = 60, cycles = 50, reads = 3, queues = 8;
  int herd_chunks = 12, herd_workers = 3, budget = 120;
  for (int i = 1; i < argc; i++) {
    std::string a = argv[i];
    auto next_val = [&]() -> int {
      CHECK(i + 1 < argc, "missing value for %s", a.c_str());
      return atoi(argv[++i]);
    };
    if      (a == "--size")         N = next_val();
    else if (a == "--cycles")       cycles = next_val();
    else if (a == "--reads")        reads = next_val();
    else if (a == "--queues")       queues = next_val();
    else if (a == "--herd-chunks")  herd_chunks = next_val();
    else if (a == "--herd-workers") herd_workers = next_val();
    else if (a == "--budget")       budget = next_val();
    else { fprintf(stderr, "unknown flag: %s\n", a.c_str()); return 64; }
  }
  g_scenario = "args";

  signal(SIGALRM, on_timeout);
  alarm(budget);                                     // hard watchdog over the whole run

  printf("buckygen-stress: size=%d cycles=%d reads=%d queues=%d herd=%dx%d budget=%ds\n",
         N, cycles, reads, queues, herd_chunks, herd_workers, budget);

  scenario_churn(N, cycles, reads);
  scenario_concurrent(N, queues);
  scenario_sigterm_handler(N, cycles < 50 ? cycles : 50);
  scenario_double_stop(20);                          // C20: a single isomer, drains instantly
  scenario_herd(N, herd_chunks, herd_workers);
  scenario_parent_death(N, budget / 4 < 5 ? 5 : budget / 4);

  alarm(0);
  printf("buckygen-stress: ALL SCENARIOS PASSED\n");
  return 0;
}
