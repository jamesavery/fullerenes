#include "libgraph/isomerdb.hh"


typedef struct {
  uint8_t RSPI[11]; 		// RSPI[2-12]
  struct {
    uint8_t Ne : 4;
    uint8_t Delta : 4;
  } HOMO;
  uint16_t iHLgap;		// 0.xxxxx -- HLgap = iHLgap/2^16
} DBEntry;

int main(int ac, char **av)
{
  if(ac<=3) return -1;
  
  int N         = strtol(av[1],0,0);
  bool IPR = strtol(av[2],0,0);
  string extension = av[3];

  IsomerDB DB(IsomerDB::readBinary(N,IPR));

  bool nontrivial = (extension == "-nontrivial");

  string filename = string(FULLERENE_DATABASE)+"/binary_new/"+(IPR?"IPR":"All")+extension+"/c"+to_string(N)+".smallbin";

  FILE *f = fopen(filename.c_str(),"wb");
  if(!f){
    perror(filename.c_str());
    exit(-1);
  }

  printf("Nisomers: %d,%d,%d\n",int(IsomerDB::number_isomers(N,nontrivial?"Nontrivial":"Any",IPR)),DB.Nisomers,int(DB.entries.size()));


  for(int i=0;i<DB.Nisomers;i++){
    IsomerDB::Entry e(DB.entries[i]);
    DBEntry be;
    be.HOMO.Ne = e.NeHOMO;
    be.HOMO.Delta = e.NedgeHOMO-e.NeHOMO/2;
    be.iHLgap = round(65536.*e.HLgap);

    for(int i=1;i<12;i++) be.RSPI[i-1] = e.RSPI[i];

    assert((e.NeHOMO&1) == 0);
    assert(be.iHLgap < 65536);

    if(nontrivial && string(e.group,3) == " C1") continue;

    size_t nwr = fwrite(&be,sizeof(DBEntry),1,f);
    assert(nwr == 1);
  }
    
  fclose(f);
  return 0;
}
