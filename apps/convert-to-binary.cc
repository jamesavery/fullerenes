#include "libgraph/isomerdb.hh"

int main(int ac, char **av)
{
  if(ac<=3) return -1;
  
  int N         = strtol(av[1],0,0);
  string outdir = string(av[2]);
  bool IPR = strtol(av[3],0,0);

  IsomerDB DB(IsomerDB::readPDB(N,IPR));

  string filename = outdir+"/c"+pad_string(to_string(N),3)+(IPR?"IPR":"all")+".bin";

  if(DB.writeBinary(filename) != true) return -2;
    
  return 0;
}
