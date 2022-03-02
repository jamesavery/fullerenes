#include <sys/stat.h>
#include <limits.h>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <thread>
#include <future>
#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/gpu/isomerspace_tutte.hh"
#include "numeric"

using namespace std;
using namespace std::chrono;

#define TEST_SAMPLES 2000

bool is_equal(const vector<coord2d>& a, const vector<coord2d> &b){
    bool result = true;
    double epsilon = 1e-4;
    result &= (a.size() == b.size());
    for (size_t i = 0; i < a.size(); i++)
    {
        result &= fabs(a[i].first - b[i].first) < epsilon;
        result &= fabs(a[i].second - b[i].second) < epsilon;
        if(!result) break;
    }
    return result;
}

int main(int ac, char **argv)
{
    if(ac<2){
    fprintf(stderr,"Syntax: %s <N:int> [output_dir] [IPR:0|1] [only_nontrivial:0|1]\n",argv[0]);
    return -1;
    }
    int N                = strtol(argv[1],0,0);     // Argument 1: Number of vertices N

    string output_dir     = ac>=3? argv[2] : "output";    // Argument 2: directory to output files to
    int IPR               = ac>=4? strtol(argv[3],0,0):0; // Argument 3: Only generate IPR fullerenes?
    int only_nontrivial   = ac>=5? strtol(argv[4],0,0):0; // Argument 4: Only generate fullerenes with nontrivial symmetry group?
    int n_best_candidates = ac>=6? strtol(argv[5],0,0):100; // Argument 5: How many best fullerne candidates do you want to store? 

    IsomerspaceTutte tutte_kernel   = IsomerspaceTutte(N);
    FullereneDual dualG;
    BuckyGen::buckygen_queue Q = BuckyGen::start(N,IPR,only_nontrivial);  
    size_t I = 0;
    bool more_to_generate = true;
    queue<FullereneGraph> fullerene_graphs;


    for (I; I < TEST_SAMPLES && more_to_generate; I++)
    {
        more_to_generate &= BuckyGen::next_fullerene(Q,dualG);
        if (!more_to_generate){break;}
        dualG.update();   		        // Update triangles
        FullereneGraph   G = dualG.dual_graph();  // Construct fullerene graph
        fullerene_graphs.push(G);
        fullerene_graphs.back().layout2d = G.tutte_layout();
        tutte_kernel.insert_isomer(G,I);
    }
    tutte_kernel.update_batch();
    while (tutte_kernel.get_batch_size()!=0 || !tutte_kernel.insert_queue.empty())
    {
        tutte_kernel.tutte_layout();
        tutte_kernel.check_batch();
        tutte_kernel.update_batch();
    }


    while (!tutte_kernel.output_queue.empty())
    {   

        if(!is_equal(fullerene_graphs.front().layout2d, tutte_kernel.output_queue.front().second.layout2d)){
            cout << fullerene_graphs.front().layout2d << std::endl;
            cout << tutte_kernel.output_queue.front().second.layout2d << std::endl;
            assert(false);
        }
        
        
        tutte_kernel.output_queue.pop(); fullerene_graphs.pop();
    }

}