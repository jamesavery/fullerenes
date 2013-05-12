#include <sqlite3.h>
#include "libgraph/isomerdb.hh"

int main(int ac, char **av)
{
  if(ac<=2) return -1;
  
  int N         = strtol(av[1],0,0);
  string outdir = string(av[2]);

  IsomerDB DB(IsomerDB::readPDB(N));

  string filename = outdir+"/c"+pad_string(to_string(N),3)+(DB.IPR?"IPR":"all")+".bin";

  if(DB.writeBinary(filename) != true) return -2;
    
  return 0;
}
