#include <iostream>
#include <cmath>
#include <cstdio>
#include <unistd.h>
#include <fcntl.h>
#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/polyhedron.hh"

using namespace std;

// Redirect stdout to /dev/null during Fortran calls (optimizer writes to unit 6 = stdout)
struct SuppressStdout {
  int saved_fd;
  int devnull;
  SuppressStdout() {
    fflush(stdout);
    saved_fd = dup(STDOUT_FILENO);
    devnull = open("/dev/null", O_WRONLY);
    dup2(devnull, STDOUT_FILENO);
  }
  ~SuppressStdout() {
    fflush(stdout);
    dup2(saved_fd, STDOUT_FILENO);
    close(saved_fd);
    close(devnull);
  }
};

int main(int ac, char **av)
{
  if(ac < 3){
    fprintf(stderr, "Usage: %s <rspi_output> <geometry_output>\n", av[0]);
    return 1;
  }

  const int N = 44;
  const char* rspi_path = av[1];
  const char* geom_path = av[2];

  FILE* f_rspi = fopen(rspi_path, "wb");
  FILE* f_geom = fopen(geom_path, "wb");
  if(!f_rspi || !f_geom){
    perror("fopen");
    return 1;
  }

  Triangulation dualG;
  BuckyGen::buckygen_queue Q = BuckyGen::start(N, false, false);

  int count = 0, written = 0;
  while(BuckyGen::next_fullerene(Q, dualG)){
    count++;
    dualG.update();

    // Extract spiral and convert to pentagon indices
    vector<int> spiral_code;
    jumplist_t jumps;
    dualG.get_spiral(spiral_code, jumps);

    uint8_t rspi[12];
    int pi = 0;
    for(int i = 0; i < (int)spiral_code.size() && pi < 12; i++)
      if(spiral_code[i] == 5)
        rspi[pi++] = (uint8_t)i;

    // Build fullerene graph and optimize geometry
    FullereneGraph G = dualG.dual_graph();
    G.layout2d = G.tutte_layout();

    vector<coord3d> points = G.zero_order_geometry();
    {
      SuppressStdout quiet;
      points = G.optimized_geometry(points);
    }

    // Check for NaN
    bool valid = true;
    for(const auto& p : points)
      if(isnan(p[0]) || isnan(p[1]) || isnan(p[2])){ valid = false; break; }
    if(!valid){
      fprintf(stderr, "Isomer %d failed optimization, skipping\n", count);
      continue;
    }

    // Center
    coord3d cm;
    for(auto& p : points) cm += p;
    cm /= (double)N;
    for(auto& p : points) p -= cm;

    // Write RSPI: 12 uint8_t
    fwrite(rspi, 1, 12, f_rspi);

    // Write geometry: 44*3 float32
    float coords[44*3];
    for(int v = 0; v < N; v++){
      coords[v*3+0] = (float)points[v][0];
      coords[v*3+1] = (float)points[v][1];
      coords[v*3+2] = (float)points[v][2];
    }
    fwrite(coords, sizeof(float), 44*3, f_geom);

    written++;
  }

  fclose(f_rspi);
  fclose(f_geom);
  fprintf(stderr, "Wrote %d of %d C%d isomers to %s and %s\n",
          written, count, N, rspi_path, geom_path);
  return 0;
}
