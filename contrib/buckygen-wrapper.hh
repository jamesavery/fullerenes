namespace BuckyGen {
  
  typedef struct {
    pid_t pid;
    int   qid;
    int   Nvertices;
    int   IPR;
  } buckygen_queue;

  buckygen_queue start(int N, int IPR);
  void stop(const buckygen_queue& Q);

  bool next_fullerene(const buckygen_queue& Q, Graph& G);
}
