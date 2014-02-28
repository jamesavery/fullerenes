#include "libgraph/planargraph.hh"
#include "libgraph/cubicgraph.hh"
#include "libgraph/fullerenegraph.hh"

using namespace std;

int main(int ac, char **av)
{
  vector<int> rspi(12);
  FullereneGraph::jumplist_t jumps;

  if(ac != 14) return -1;

  int N = atol(av[1]);
  cout << N << endl;
  for(int i=0;i<12;i++) rspi[i] = strtol(av[2+i],0,0)-1;

  cout << "Attempting to create graph from spiral indices " << rspi << endl;
 
  FullereneGraph g(N, rspi, jumps);
  cout << "done" << endl;

  vector<int> pent_dist_vec = g.pentagon_distance_mtx();

  cout << "{";
  for(int i=0; i<12; ++i)
//     cout << vector<int>(&pent_dist_mtx[12*i],&pent_dist_mtx[13*i]) << (i+1<12?",\n":"");
     cout << vector<int>(&pent_dist_vec[12*i],&pent_dist_vec[12*(i+1)]) << (i+1<12?",\n ":"");
  cout << "}" << endl;

// store mtx in 12x12 array
  int pd[12][12];
  for(int i=0; i<12; ++i){
    copy(pent_dist_vec.begin()+12*i, pent_dist_vec.begin()+12*(i+1), pd[i]);
  }
  
//  for(int i=0; i<12; ++i){
//    for(int j=0; j<12; ++j){
//      cout << pd[i][j] << ", ";
//    }
//    cout << endl;
//  }

  // find C_6 axes
  int c6_count = 0;
  
  // find all ways to choose 6 distinct numbers out of 12
  // this is one set of six vertices that is subjected to cyclic shift
  // fix 0 in the first position of the first ring to avoid double counting
  int i = 0;
  for (int j=1; j<12; ++j){
    if(i==j) continue;
    for(int k=1; k<12; ++k){
      if(k==i || k==j) continue;
      for(int l=1; l<12; ++l){
        if(l==i || l==j || l==k) continue;
        for(int m=1; m<12; ++m){
          if(m==i || m==j || m==k || m==l) continue;
          for(int n=1; n<12; ++n){
            if(n==i || n==j || n==k || n==l || n==m) continue;

            // find all permutations of the remaining 6 numbers
            // this is the second set of six vertices that is subjected to cyclic shift
            int foo[]= {0,1,2,3,4,5,6,7,8,9,10,11};
            vector<int> o (foo, foo+12);
            o.erase(std::remove(o.begin(), o.end(), i), o.end());
            o.erase(std::remove(o.begin(), o.end(), j), o.end());
            o.erase(std::remove(o.begin(), o.end(), k), o.end());
            o.erase(std::remove(o.begin(), o.end(), l), o.end());
            o.erase(std::remove(o.begin(), o.end(), m), o.end());
            o.erase(std::remove(o.begin(), o.end(), n), o.end());
//        cout << i << ", " << j << ", " << k << ", " << l << ", " << m << ", " << n << endl;
     
            // distances between pairs of vertices in the first set, before and after the shift
            if (!(pd[i][i]==pd[n][n] && pd[i][j]==pd[n][i] && pd[i][k]==pd[n][j] && pd[i][l]==pd[n][k] && pd[i][m]==pd[n][l] && pd[i][n]==pd[n][m])) continue;
            if (!(pd[j][i]==pd[i][n] && pd[j][j]==pd[i][i] && pd[j][k]==pd[i][j] && pd[j][l]==pd[i][k] && pd[j][m]==pd[i][l] && pd[j][n]==pd[i][m])) continue;
            if (!(pd[k][i]==pd[j][n] && pd[k][j]==pd[j][i] && pd[k][k]==pd[j][j] && pd[k][l]==pd[j][k] && pd[k][m]==pd[j][l] && pd[k][n]==pd[j][m])) continue;
            if (!(pd[l][i]==pd[k][n] && pd[l][j]==pd[k][i] && pd[l][k]==pd[k][j] && pd[l][l]==pd[k][k] && pd[l][m]==pd[k][l] && pd[l][n]==pd[k][m])) continue;
            if (!(pd[m][i]==pd[l][n] && pd[m][j]==pd[l][i] && pd[m][k]==pd[l][j] && pd[m][l]==pd[l][k] && pd[m][m]==pd[l][l] && pd[m][n]==pd[l][m])) continue;
            if (!(pd[n][i]==pd[m][n] && pd[n][j]==pd[m][i] && pd[n][k]==pd[m][j] && pd[n][l]==pd[m][k] && pd[n][m]==pd[m][l] && pd[n][n]==pd[m][m])) continue;
            do {
              if(
                  // distances between pairs of vertices in the first and second set, before and after the shift
                  pd[i][o[0]]==pd[n][o[5]] && pd[i][o[1]]==pd[n][o[0]] && pd[i][o[2]]==pd[n][o[1]] && pd[i][o[3]]==pd[n][o[2]] && pd[i][o[4]]==pd[n][o[3]] && pd[i][o[5]]==pd[n][o[4]] &&
                  pd[j][o[0]]==pd[i][o[5]] && pd[j][o[1]]==pd[i][o[0]] && pd[j][o[2]]==pd[i][o[1]] && pd[j][o[3]]==pd[i][o[2]] && pd[j][o[4]]==pd[i][o[3]] && pd[j][o[5]]==pd[i][o[4]] &&  
                  pd[k][o[0]]==pd[j][o[5]] && pd[k][o[1]]==pd[j][o[0]] && pd[k][o[2]]==pd[j][o[1]] && pd[k][o[3]]==pd[j][o[2]] && pd[k][o[4]]==pd[j][o[3]] && pd[k][o[5]]==pd[j][o[4]] &&  
                  pd[l][o[0]]==pd[k][o[5]] && pd[l][o[1]]==pd[k][o[0]] && pd[l][o[2]]==pd[k][o[1]] && pd[l][o[3]]==pd[k][o[2]] && pd[l][o[4]]==pd[k][o[3]] && pd[l][o[5]]==pd[k][o[4]] &&  
                  pd[m][o[0]]==pd[l][o[5]] && pd[m][o[1]]==pd[l][o[0]] && pd[m][o[2]]==pd[l][o[1]] && pd[m][o[3]]==pd[l][o[2]] && pd[m][o[4]]==pd[l][o[3]] && pd[m][o[5]]==pd[l][o[4]] &&  
                  pd[n][o[0]]==pd[m][o[5]] && pd[n][o[1]]==pd[m][o[0]] && pd[n][o[2]]==pd[m][o[1]] && pd[n][o[3]]==pd[m][o[2]] && pd[n][o[4]]==pd[m][o[3]] && pd[n][o[5]]==pd[m][o[4]] &&
                  // distances between pairs of vertices in the second set, before and after the shift
                  pd[o[0]][o[0]]==pd[o[5]][o[5]] && pd[o[0]][o[1]]==pd[o[5]][o[0]] && pd[o[0]][o[2]]==pd[o[5]][o[1]] && pd[o[0]][o[3]]==pd[o[5]][o[2]] && pd[o[0]][o[4]]==pd[o[5]][o[3]] && pd[o[0]][o[5]]==pd[o[5]][o[4]] &&  
                  pd[o[1]][o[0]]==pd[o[0]][o[5]] && pd[o[1]][o[1]]==pd[o[0]][o[0]] && pd[o[1]][o[2]]==pd[o[0]][o[1]] && pd[o[1]][o[3]]==pd[o[0]][o[2]] && pd[o[1]][o[4]]==pd[o[0]][o[3]] && pd[o[1]][o[5]]==pd[o[0]][o[4]] &&  
                  pd[o[2]][o[0]]==pd[o[1]][o[5]] && pd[o[2]][o[1]]==pd[o[1]][o[0]] && pd[o[2]][o[2]]==pd[o[1]][o[1]] && pd[o[2]][o[3]]==pd[o[1]][o[2]] && pd[o[2]][o[4]]==pd[o[1]][o[3]] && pd[o[2]][o[5]]==pd[o[1]][o[4]] &&  
                  pd[o[3]][o[0]]==pd[o[2]][o[5]] && pd[o[3]][o[1]]==pd[o[2]][o[0]] && pd[o[3]][o[2]]==pd[o[2]][o[1]] && pd[o[3]][o[3]]==pd[o[2]][o[2]] && pd[o[3]][o[4]]==pd[o[2]][o[3]] && pd[o[3]][o[5]]==pd[o[2]][o[4]] &&  
                  pd[o[4]][o[0]]==pd[o[3]][o[5]] && pd[o[4]][o[1]]==pd[o[3]][o[0]] && pd[o[4]][o[2]]==pd[o[3]][o[1]] && pd[o[4]][o[3]]==pd[o[3]][o[2]] && pd[o[4]][o[4]]==pd[o[3]][o[3]] && pd[o[4]][o[5]]==pd[o[3]][o[4]] &&  
                  pd[o[5]][o[0]]==pd[o[4]][o[5]] && pd[o[5]][o[1]]==pd[o[4]][o[0]] && pd[o[5]][o[2]]==pd[o[4]][o[1]] && pd[o[5]][o[3]]==pd[o[4]][o[2]] && pd[o[5]][o[4]]==pd[o[4]][o[3]] && pd[o[5]][o[5]]==pd[o[4]][o[4]]){
                 // increment if distance between two vertices is unaffected by mapping of graph on itself
                 ++c6_count;
                 cout << i << ", " << j << ", " << k << ", " << l << ", " << m << ", " << n << o << endl;
               }
            } while ( std::next_permutation(o.begin()+1, o.end()) );//leave the first in the first position to avoid double counting
          }
        }
      }
    }
  }
  cout << c6_count << " C6 elements found " << endl;

  return 0;
}

