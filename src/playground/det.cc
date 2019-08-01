#include <stdio.h>
#include <iostream>
#include <vector>
#include "fullerenes/fullerenegraph.hh"
#include "fullerenes/isomerdb.hh"

using namespace std;
string pad_string(const string& s, int length, char padchar)
{
  string result(s);
  int pad = length - s.size();
  string padstring;
  for(int i=0;i<pad;i++) padstring += padchar;
  return padstring+result;
}

int64_t det(vector< vector<int> > A, int n);
double lu_det(const vector<double> &A, int N);
void perfmatch_dfs(map<dedge_t,int>& faceEdge, const vector<face_t>& faces, 
		   map<dedge_t,int>& matrix, vector<bool>& faceSum, vector<bool>& visited, const dedge_t& e);


int main(int ac, char **av)
{
  if(ac<2) return -1;
  int N = strtol(av[1],0,0);

  IsomerDB DB(IsomerDB::readBinary("database/binary/c"+pad_string(to_string(N),3)+"all.bin"));
 
  //  printf("Number of perfect matchings for all %d isomers of C%d\n",DB.Nisomers,N);
  for(int isomer=0;isomer<DB.Nisomers;isomer++){
    IsomerDB::Entry e(DB.entries[isomer]);
    vector<int> rspi(e.RSPI,e.RSPI+12);
    for(int i=0;i<12;i++) rspi[i]--;

    FullereneGraph G(N,rspi);
    G.layout2d = G.tutte_layout();

    cout << G.count_perfect_matchings() << endl;
  } 
  return 0;
}
