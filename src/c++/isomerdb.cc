#include "isomerdb.hh"
// TODO: Hov! Isomer count is wrong when reading text database. Find and fix!

// TODO: Make C++98 compatible to support old compilers or make optional

string IsomerDB::database_path = FULLERENE_DATABASE;

// Metadata for Cn isomers
vector<size_t> IsomerDB::Nisomers_data[2]                 = {
  {1,0,1,1,2,3,6,6,15,17,40,45,89,116,199,271,437,580,924,1205,1812,2385,3465,4478,6332,8149,11190,14246,19151,24109,31924,39718,51592,63761,81738,99918,126409,153493,191839,231017,285914,341658,419013,497529,604217,713319,860161,1008444,1207119,1408553,1674171,1942929,2295721,2650866,3114236,3580637,4182071,4787715,5566949,6344698,7341204,8339033,9604411,10867631,12469092,14059174,16066025,18060979,20558767,23037594,26142839,29202543,33022573,36798433,41478344,46088157,51809031,57417264,64353269,71163452,79538751,87738311,97841183,107679717,119761075,131561744,145976674,159999462,177175687,193814658,214127742},
  {1,0,0,0,0,1,1,1,2,5,7,9,24,19,35,46,86,134,187,259,450,616,823,1233,1799,2355,3342,4468,6063,8148,10774,13977,18769,23589,30683,39393,49878,62372,79362,98541,121354,151201,186611,225245,277930,335569,404667,489646,586264,697720,836497,989495,1170157,1382953,1628029,1902265,2234133,2601868,3024383,3516365,4071832,4690880,5424777,6229550,7144091,8187581,9364975,10659863,12163298,13809901,15655672}
};

vector< vector<string> > IsomerDB::symmetries_data[2]     = {{{" Ih"},{},{"D6d"},{"D3h"},{" D2"," Td"},{"C2v","D5h"},{" C2"," D2"," D3","D3d","D3h"},{" C2"," Cs","C3v"},{" C1"," C2"," Cs"," D2","C2v","D2d","D3h","D6h"},{" C1"," C2"," D3","C2v","C3v","D3h"},{" C1"," C2"," C3"," Cs"," D2"," Td","C2v","C3v","D2h","D5d"},{" C1"," C2"," Cs"," D3","C2v"},{"  T"," C1"," C2"," Cs"," D2"," D3"," S4","C2v","D3d","D3h"},{" C1"," C2"," C3"," Cs","C2v"},{" C1"," C2"," Cs"," D2"," D3","C2h","C2v","D2h","D6d"},{" C1"," C2"," C3"," Cs"," D3","C2v","C3v","D3h","D5h"},{"  T"," C1"," C2"," C3"," Cs"," D2","C2h","C2v","C3v","D2d","D2h"},{" C1"," C2"," Cs"," D3","C2v","D3h"},{" C1"," C2"," C3"," Ci"," Cs"," D2"," D3"," Td","C2h","C2v","D2d","D2h","D3d"},{" C1"," C2"," C3"," Cs","C2v","C3v"},{" C1"," C2"," Cs"," D2"," D3"," D5"," Ih"," S4","C2h","C2v","C3v","D2d","D2h","D5d","D6h"},{" C1"," C2"," C3"," Cs"," D3","C2v","C3h","C3v","D3h"},{" C1"," C2"," C3"," Ci"," Cs"," D2","C2h","C2v","C3v","D2d","D2h"},{" C1"," C2"," Cs"," D3","C2v","C3v"},{"  T"," C1"," C2"," C3"," Cs"," D2"," D3"," S4"," S6"," Td","C2h","C2v","C3h","D2d","D2h","D3d","D3h"},{" C1"," C2"," C3"," Cs","C2v","C3v","D5h"},{" C1"," C2"," C3"," Ci"," Cs"," D2"," D3"," D6","C2h","C2v","C3v","D2d","D2h","D6d"},{" C1"," C2"," C3"," Cs"," D3","C2v","C3h","C3v","D3h"},{"  T"," C1"," C2"," C3"," Ci"," Cs"," D2"," S4"," Td","C2h","C2v","C3v","D2d"},{" C1"," C2"," C3"," Cs"," D3","C2v","D3h"},{" C1"," C2"," C3"," Ci"," Cs"," D2"," D3"," Ih"," S4","C2h","C2v","C3v","D3d","D3h","D5d","D5h"},{" C1"," C2"," C3"," Cs","C2v","C3v"},{" C1"," C2"," C3"," Ci"," Cs"," D2"," D3"," S4"," Td","C2h","C2v","C3v","D2d","D2h","D3d","D3h","D6h"},{" C1"," C2"," C3"," Cs"," D3","C2v","C3v","D3h"},{"  T"," C1"," C2"," C3"," Ci"," Cs"," D2"," S4","C2h","C2v","C3v","D2d","D2h"},{" C1"," C2"," C3"," Cs"," D3"," D5","C2v","D3h","D5h"},{"  T"," C1"," C2"," C3"," Ci"," Cs"," D2"," D3"," S4"," Td"," Th","C2h","C2v","C3h","C3v","D2d","D2h","D3d","D3h"},{" C1"," C2"," C3"," Cs","C2v","C3v"},{" C1"," C2"," C3"," Ci"," Cs"," D2"," D3","C2h","C2v","C3v","D2d","D2h","D3d","D3h","D6d","D6h"},{" C1"," C2"," C3"," Cs"," D3","C2v","C3h","C3v","D3h"},{"  T"," C1"," C2"," C3"," Ci"," Cs"," D2"," D5"," S4"," Td","C2h","C2v","C3v","D2d","D2h","D5d"},{" C1"," C2"," C3"," Cs"," D3","C2v","C3v"},{"  T"," C1"," C2"," C3"," Ci"," Cs"," D2"," D3"," S4","C2h","C2v","C3v","D2d","D2h","D3d","D3h"},{" C1"," C2"," C3"," Cs","C2v","C3v"},{" C1"," C2"," C3"," Ci"," Cs"," D2"," D3"," D6"," S4","C2h","C2v","C3h","C3v","D2d","D2h","D3d","D3h","D6h"},{" C1"," C2"," C3"," Cs"," D3"," D5","C2v","C3v","D3h","D5h"},{" C1"," C2"," C3"," Ci"," Cs"," D2"," S4"," Td","C2h","C2v","C3v","D2d","D2h"},{" C1"," C2"," C3"," Cs"," D3","C2v","C3v","D3h"},{"  T"," C1"," C2"," C3"," Ci"," Cs"," D2"," D3"," S4"," S6"," Th","C2h","C2v","C3h","C3v","D2d","D2h","D3d","D3h"},{" C1"," C2"," C3"," Cs","C2v","C3v"},{" C1"," C2"," C3"," Ci"," Cs"," D2"," D3"," D5"," D6"," S4"," Td","C2h","C2v","C3v","D2d","D2h","D5d","D5h","D6d"}},

 {{" Ih"},{},{},{},{},{"D5h"},{"D6d"},{"D3h"},{" D2"," Td"},{" D3","C2v","D3h"},{" D2"," D3"," Ih","C2v","D5d","D5h"},{" C2"," Cs","C2v","C3v"},{" C1"," C2"," Cs"," D2"," Td","C2v","D2d","D3d","D6h"},{" C1"," C2"," C3"," Cs"," D3","C2v"},{"  T"," C1"," C2"," Cs"," D2","C2v"},{" C1"," C2"," Cs","C2v","D5h"},{"  T"," C1"," C2"," C3"," Cs"," D2"," D3","C2v","D2h"},{" C1"," C2"," C3"," Cs","C2v","C3v"},{" C1"," C2"," Cs"," D2"," D3","C2v","C3v","D2d","D2h","D3d","D3h","D6d","D6h"},{" C1"," C2"," C3"," Cs"," D3","C2v"},{"  T"," C1"," C2"," C3"," Cs"," D2"," D5","C2v","D2d","D5d"},{" C1"," C2"," Cs"," D3","C2v","C3v"},{" C1"," C2"," C3"," Cs"," D2"," D3","C2v","D2d","D2h","D3d","D3h"},{" C1"," C2"," C3"," Cs","C2v","C3v"},{" C1"," C2"," Cs"," D2"," D3"," S4","C2h","C2v","C3v","D2d","D2h","D3d","D3h","D6h"},{" C1"," C2"," C3"," Cs"," D3"," D5","C2v","D5h"},{" C1"," C2"," C3"," Cs"," D2"," Td","C2h","C2v","C3v","D2h"},{" C1"," C2"," C3"," Cs"," D3","C2v","C3v","D3h"},{"  T"," C1"," C2"," C3"," Cs"," D2"," D3"," Th","C2h","C2v","C3h","D2d","D2h"},{" C1"," C2"," C3"," Cs","C2v","C3v"},{" C1"," C2"," C3"," Ci"," Cs"," D2"," D3"," D5"," D6"," Td","C2h","C2v","C3v","D2d","D2h","D5d","D5h","D6d"},{" C1"," C2"," C3"," Cs"," D3","C2v"},{"  T"," C1"," C2"," C3"," Ci"," Cs"," D2"," S4","C2h","C2v","C3v","D2d","D2h"},{" C1"," C2"," C3"," Cs"," D3","C2v"},{" C1"," C2"," C3"," Ci"," Cs"," D2"," D3"," S6","C2h","C2v","D2d","D3d","D3h"},{" C1"," C2"," C3"," Cs"," D5","C2v","C3v","D5h"},{"  T"," C1"," C2"," C3"," Cs"," D2"," D3"," D6"," S4","C2h","C2v","C3v","D2d","D2h","D3d","D3h","D6h"},{" C1"," C2"," C3"," Cs"," D3","C2v","C3h","C3v","D3h"},{"  T"," C1"," C2"," C3"," Ci"," Cs"," D2"," S4","C2h","C2v","C3v","D2d","D2h"},{" C1"," C2"," C3"," Cs"," D3","C2v","C3v"},{"  I","  T"," C1"," C2"," C3"," Ci"," Cs"," D2"," D3"," D5"," S4","C2h","C2v","C3h","D5d"},{" C1"," C2"," C3"," Cs","C2v","C3v"},{" C1"," C2"," C3"," Ci"," Cs"," D2"," D3"," D6"," S4","C2h","C2v","C3v","D2d","D2h","D3d","D3h","D6d","D6h"},{" C1"," C2"," C3"," Cs"," D3","C2v","D3h"},{"  T"," C1"," C2"," C3"," Ci"," Cs"," D2"," S4","C2h","C2v","C3v","D2d","D2h"},{" C1"," C2"," C3"," Cs"," D3"," D5","C2v","C3v","D3h","D5h"},{"  T"," C1"," C2"," C3"," Ci"," Cs"," D2"," D3"," Th","C2h","C2v","C3v","D2d","D3d","D3h"},{" C1"," C2"," C3"," Cs","C2v","C3v"},{"  T"," C1"," C2"," C3"," Ci"," Cs"," D2"," D3"," D6"," S4","C2h","C2v","C3v","D2d","D2h","D3d","D3h","D6h"},{" C1"," C2"," C3"," Cs"," D3","C2v","C3h"},{" C1"," C2"," C3"," Ci"," Cs"," D2"," D5"," S4"," Td","C2h","C2v","C3v","D2d","D2h","D5d","D5h"},{" C1"," C2"," C3"," Cs"," D3","C2v","C3v","D3h"},{"  T"," C1"," C2"," C3"," Ci"," Cs"," D2"," D3"," S4"," S6"," Td","C2h","C2v","C3h","D2d","D2h","D3d"},{" C1"," C2"," C3"," Cs","C2v","C3v"},{" C1"," C2"," C3"," Ci"," Cs"," D2"," D3"," D6"," S4"," Td","C2h","C2v","C3v","D2d","D2h","D3d","D6d"},{" C1"," C2"," C3"," Cs"," D3"," D5","C2v","C3v","D3h","D5h"},{"  T"," C1"," C2"," C3"," Ci"," Cs"," D2"," S4"," Td","C2h","C2v","D2d","D2h"},{" C1"," C2"," C3"," Cs"," D3","C2v","C3v","D3h"},{"  T"," C1"," C2"," C3"," Ci"," Cs"," D2"," D3"," S4","C2h","C2v","C3h","C3v","D2d","D2h","D3d","D3h"},{" C1"," C2"," C3"," Cs","C2v","C3v"},{" C1"," C2"," C3"," Ci"," Cs"," D2"," D3"," D5"," D6"," Ih"," S4","C2h","C2v","C3v","D2d","D2h","D3d","D3h","D5d","D6h"},{" C1"," C2"," C3"," Cs"," D3","C2v","C3h"},{" C1"," C2"," C3"," Ci"," Cs"," D2"," S4","C2h","C2v","C3v","D2d","D2h"},{" C1"," C2"," C3"," Cs"," D3","C2v","C3h","C3v","D3h"},{"  T"," C1"," C2"," C3"," Ci"," Cs"," D2"," D3"," S4"," S6"," Th","C2h","C2v","C3v","D2d","D2h","D3d","D3h"},{" C1"," C2"," C3"," Cs"," D5","C2v","C3v","D5h"},{" C1"," C2"," C3"," Ci"," Cs"," D2"," D3"," D6"," S4","C2h","C2v","C3v","D2d","D2h","D3d","D6d","D6h"},{" C1"," C2"," C3"," Cs"," D3","C2v","C3h","C3v","D3h"},{"  T"," C1"," C2"," C3"," Ci"," Cs"," D2"," S4","C2h","C2v","C3v","D2d","D2h"},{" C1"," C2"," C3"," Cs"," D3","C2v","C3v","D3h"},{"  T"," C1"," C2"," C3"," Ci"," Cs"," D2"," D3"," D5"," S4"," Th","C2h","C2v","C3v","D2d","D2h","D3d","D3h","D5d","D5h"}}
};

vector< vector<size_t> > IsomerDB::symmetry_count_data[2] = {{{1},{},{1},{1},{1,1},{2,1},{2,1,1,1,1},{3,2,1},{2,4,2,2,1,2,1,1},{7,5,1,2,1,1},{8,14,1,7,3,1,2,1,1,2},{23,11,6,1,4},{1,42,22,7,6,2,1,3,3,2},{69,22,2,19,4},{117,52,16,5,1,1,3,2,2},{195,37,2,25,2,6,1,1,2},{1,307,78,3,26,9,1,3,2,5,2},{470,62,38,1,8,1},{700,135,3,1,49,10,6,1,2,13,1,1,2},{1037,98,4,58,6,2},{1508,189,67,19,3,1,1,2,4,9,1,4,1,1,2},{2135,142,4,80,4,16,1,1,2},{2990,316,8,2,118,17,4,4,4,1,1},{4134,211,112,2,18,1},{1,5714,411,5,122,28,10,2,1,1,7,21,1,2,2,3,1},{7634,300,8,186,14,5,2},{10304,619,1,3,190,24,3,1,7,26,1,4,5,2},{13557,414,9,237,6,18,1,1,3},{1,18005,800,14,2,246,45,4,1,11,14,5,3},{23197,557,2,312,3,35,3},{30280,1146,15,5,371,39,12,1,1,16,28,2,2,2,2,2},{38548,742,15,380,28,5},{49590,1436,1,9,434,59,6,3,1,4,29,1,9,6,1,1,2},{62212,976,15,505,9,36,5,3},{1,79033,1945,24,10,596,52,1,16,43,8,4,5},{97936,1266,3,655,4,1,50,1,2},{1,123141,2412,20,12,646,80,19,5,1,1,13,38,2,2,4,4,5,3},{150939,1603,26,879,42,4},{187505,3200,4,20,972,70,9,16,28,2,3,3,1,1,3,2},{227934,2029,27,952,13,58,1,1,2},{2,280730,3801,40,14,1093,114,2,9,1,28,66,5,5,2,2},{337808,2542,3,1228,6,67,4},{1,412339,4954,37,29,1413,89,24,1,23,82,2,7,7,4,1},{492768,3126,38,1541,48,8},{596532,5872,5,26,1501,145,10,1,9,29,65,1,1,10,3,2,3,2},{707444,3845,42,1872,14,1,90,5,4,2},{850295,7403,59,41,2147,116,1,1,28,54,8,2,6},{1001569,4684,7,2079,9,88,5,3},{2,1195728,8713,49,44,2238,169,34,7,1,1,22,76,4,2,12,10,5,2},{1400184,5610,54,2609,84,12},{1660007,10787,10,52,2921,157,11,2,2,3,1,51,139,5,8,8,4,1,2}},

{{1},{},{},{},{},{1},{1},{1},{1,1},{1,2,2},{1,1,1,2,1,1},{3,3,1,2},{1,5,5,4,1,4,2,1,1},{6,6,1,3,1,2},{1,11,7,11,2,3},{16,16,6,7,1},{1,38,26,1,8,4,5,2,1},{89,26,3,13,2,1},{108,43,14,8,3,3,1,1,1,1,1,2,1},{169,49,3,30,3,5},{1,336,62,3,31,9,1,5,1,1},{488,73,44,2,8,1},{644,123,1,26,11,4,8,2,2,1,1},{1054,105,9,54,8,3},{1479,201,72,21,2,1,3,9,1,5,1,1,2,1},{2111,168,4,56,6,1,8,1},{2950,250,6,105,20,1,2,5,2,1},{4089,258,1,94,4,18,3,1},{1,5508,386,7,112,22,10,1,2,7,2,4,1},{7670,303,11,148,13,3},{9904,612,1,1,180,29,4,1,1,1,5,25,1,1,3,3,1,1},{13295,472,12,178,7,13},{2,17751,735,15,1,200,41,3,5,10,2,3,1},{22686,625,2,242,5,29},{29275,1052,12,3,269,36,11,1,4,13,1,2,4},{38285,773,24,274,1,29,6,1},{1,48013,1369,2,346,71,8,1,6,10,35,2,3,2,5,3,1},{60767,1117,13,433,13,26,1,1,1},{1,77072,1740,23,4,408,60,2,12,32,2,3,3},{96598,1382,6,517,5,31,2},{1,2,118321,2335,24,8,535,71,17,1,4,7,26,1,1},{148844,1676,31,612,35,3},{182804,2990,4,3,665,80,8,1,1,9,27,1,5,7,1,1,3,1},{222281,2171,33,710,12,37,1},{1,273339,3581,49,11,790,98,1,9,34,3,10,4},{331710,2792,7,982,10,2,57,4,2,3},{1,398804,4681,26,15,964,91,24,1,15,39,1,2,2,1},{485184,3275,44,1099,37,7},{1,579044,5752,8,17,1203,145,10,1,5,15,41,3,10,5,1,2,1},{692268,4051,48,1280,20,51,2},{828028,6732,59,19,1461,127,1,1,1,15,39,5,3,2,3,1},{982795,5083,7,1520,10,74,3,3},{1,1159784,8415,44,21,1561,176,34,11,2,1,26,70,2,3,4,2},{1375066,5852,68,1901,58,8},{1615320,10304,15,33,2039,163,19,2,2,1,20,89,5,6,6,4,1},{1893005,7049,62,2065,19,1,60,1,2,1},{3,2219523,11833,79,30,2374,207,9,1,29,41,3,1},{2590514,8685,15,2544,12,93,3,2},{1,3006861,14377,77,41,2704,182,36,1,20,61,1,1,6,7,3,4},{3503608,9817,104,2736,87,13},{4051036,17134,12,50,3115,286,15,5,2,1,12,37,96,4,12,5,2,4,1,3},{4675522,11800,73,3357,32,93,3},{5359303,18401,90,55,3637,179,3,27,90,6,4,5},{6211666,13865,26,3874,16,94,1,5,3},{2,7116261,23225,100,64,3946,302,50,11,1,1,35,80,1,7,3,1,1},{8167004,15775,119,4567,2,110,3,1},{9331690,27698,28,91,5028,300,20,1,2,32,61,6,6,7,1,3,1},{10636469,18521,118,4586,30,133,1,2,3},{4,12126580,30655,150,78,5232,383,16,47,138,7,6,2},{13781897,21759,22,6068,18,132,2,3},{1,15612518,36317,126,97,5982,330,51,3,8,1,58,163,2,3,4,2,1,3,2}}
};



// ------------------------------ Functions ------------------------------
bool IsomerDB::writeBinary(const string filename) const
{
  FILE *f = fopen(filename.c_str(),"wb");
  if(!f){
    cerr << "Couldn't open database file " << filename << " for writing: " << strerror(errno) << ".\n";
    return false;
  }
  u_int16_t header = N | (IPR<<8) | (with_ncycham<<9);
  fwrite(&header,2,1,f);

  fwrite(&Nisomers,4,1,f);
  fwrite(&entries[0],sizeof(Entry),Nisomers,f);

  return true;
}


bool IsomerDB::writeCSV(const string filename) const
  {
    ofstream f(filename.c_str());
    if(!f){
      cerr << "Couldn't open database file " << filename << " for writing: " << strerror(errno) << ".\n";
      return false;
    }
    
    f << "\"SYMGROUP\", \"RSPI\",\"PNI\",\"HNI\",\"NEHOMO\",\"NEDGEHOMO\",\"HLGAP\",\"NCYCHAM\",\"INMR\"\r\n";
    for(int i=0;i<Nisomers;i++){
      Entry e(entries[i]);
      f << "\""<<string(e.group,3) << "\",\"" << CSVarray(e.RSPI,12) << "\",\"" << CSVarray(e.PNI,5) << "\",\"" << CSVarray(e.HNI,6) << "\","
	<< int(e.NeHOMO) << "," << int(e.NedgeHOMO) << "," << e.HLgap << "," << e.ncycham << ",\"" << CSVarray(e.INMR,6) << "\"\r\n";
    }
    return true;
  }

IsomerDB IsomerDB::readBinary(const string filename){
  FILE *f = fopen(filename.c_str(),"rb");
  if(!f){
    cerr << "Couldn't open database file " << filename << " for reading: " << strerror(errno) << ".\n";
    return IsomerDB(-1);
  }
  u_int16_t header;
  size_t nread1 = fread(&header,2,1,f);
    
  IsomerDB DB(header & 0xff, header >> 8 & 1, header >> 9 & 1);
  size_t nread2 = fread(&DB.Nisomers,4,1,f);
  DB.entries.resize(DB.Nisomers);
  size_t nread3 = fread(&DB.entries[0],sizeof(Entry),DB.Nisomers,f);

  if(!nread1 || !nread2 || !nread3) return IsomerDB(-1);
  return DB;
}
 
IsomerDB::Entry IsomerDB::getIsomer(int N, int isomer, bool IPR){
  string filename = database_path+"/binary/c"+pad_string(to_string(N),3)+(IPR?"IPR":"all")+".bin";
  FILE *f = fopen(filename.c_str(),"rb");
  if(!f){
    cerr << "Couldn't open database file " << filename << " for reading: " << strerror(errno) << ".\n";
    abort();
  }
  u_int16_t header;
  u_int32_t Nisomers;
  Entry e;

  if(fread(&header,2,1,f) != 1){ cerr << "Read error from database file " << filename << ": " << strerror(errno) << ".\n"; }
  if(fread(&Nisomers,4,1,f) != 1){ cerr << "Read error from database file " << filename << ": " << strerror(errno) << ".\n"; }
  assert(isomer <= Nisomers);
  fseek(f,(isomer-1)*sizeof(Entry),SEEK_CUR);
    
  if(fread(&e,sizeof(Entry),1,f) != 1){ cerr << "Read error from database file " << filename << ": " << strerror(errno) << ".\n"; }

  return e;
}

FullereneGraph IsomerDB::makeIsomer(int N, const Entry& e)
{
  vector<int> RSPI(e.RSPI,e.RSPI+12);
  for(int i=0;i<12;i++) RSPI[i]--;
  //    cout << "creating C"<<N<< " from spiral indices " << RSPI << endl;
  return FullereneGraph(N,RSPI);
}

IsomerDB IsomerDB::readBinary(int N, bool IPR, string extension) {
  string filename = database_path+"/binary/c"+pad_string(to_string(N),3,'0')+(IPR?"IPR":"all")+extension+string(".bin");
  return readBinary(filename);
}

IsomerDB IsomerDB::readPDB(int N, bool IPR, string extension) {
  string filename;
  if(IPR)
    filename= database_path+"/IPR/c"+pad_string(to_string(N),3,'0')+"IPR"+extension+".database";      
  else 
    filename= database_path+"/All/c"+pad_string(to_string(N),3,'0')+"all"+extension+".database";

  ifstream dbfile(filename.c_str());
  if(!dbfile){
    cerr << "Couldn't open database file " << filename << " for reading. (error: " << strerror(errno) << ")\n";
    return IsomerDB(-1);
  }
    
  string line;
  int Nread,IP,IH, pos=0;
  getline(dbfile,line);
    
  Nread = fortran_readI(line,pos,3);
  IP = fortran_readI(line,pos,1);
  IH = fortran_readI(line,pos,1);

  IsomerDB DB(Nread,IP,IH);
    
  /* 
     IP=0,IH=1:  1004 Format(A3,12I3,5I2,6I2,I2,I1,F7.5,I7,6I3)
     IP=0,IH=0:  1007 Format(A3,12I3,5I2,6I2,I2,I1,F7.5,6I3)
     IP=1,IH=1:  1008 Format(A3,12I3,3I2,I2,I1,F7.5,I7,6I3)
     IP=1,IH=0:  1009 Format(A3,12I3,3I2,I2,I1,F7.5,6I3)
  */
  while(getline(dbfile,line)){
    int pos = 0;
    Entry e;
    fortran_readA(e.group,line,pos,3);
    for(int i=0;i<12;i++) e.RSPI[i] = fortran_readI(line,pos,3);
    if(IP==0){
      for(int i=0;i<5;i++) e.PNI[i] = fortran_readI(line,pos,2);
      for(int i=0;i<3;i++) e.HNI[i] = 0;
      for(int i=3;i<6;i++) e.HNI[i] = fortran_readI(line,pos,2);
      // for(int i=0;i<6;i++) e.HNI[i] = fortran_readI(line,pos,2); // Discrepancy between util.f and db files
    } else {
      for(int i=0;i<5;i++) e.PNI[i] = 0;
      for(int i=0;i<3;i++) e.HNI[i] = 0;
      for(int i=3;i<6;i++) e.HNI[i] = fortran_readI(line,pos,2);
    }
    e.NeHOMO = fortran_readI(line,pos,2);
    e.NedgeHOMO = fortran_readI(line,pos,1); 
    e.HLgap = fortran_readF(line,pos,10);
    if(IH==1) e.ncycham = fortran_readI(line,pos,7);
    for(int i=0;i<6;i++) e.INMR[i] = fortran_readI(line,pos,3);
    DB.entries.push_back(e);
  }
    
  DB.Nisomers = DB.entries.size();
  for(int i=0;i<DB.Nisomers;i++){
    const Entry &e(DB.entries[i]);
    DB.RSPIindex[vector<int>(e.RSPI,e.RSPI+12)] = i;
  }

  return DB;
}

int64_t IsomerDB::number_isomers(int N, const string& sym, bool IPR){ 
  int Nindex = (N-(IPR?60:20))/2;
  // fprintf(stderr,"number_isomers(%d,%s,%d) : %d = %s\n",N,sym.c_str(),IPR,Nindex,
  // 	  (Nindex >= Nisomers_data[IPR].size()? "Out of range" : to_string(Nisomers_data[IPR][Nindex]).c_str())
  // 	  );
  if(Nindex >= Nisomers_data[IPR].size()) return 0;

  if(sym == "Any" || sym == "") return Nisomers_data[IPR][Nindex]; 

  if(sym == "Nontrivial"){
    size_t sum = 0;
    if(Nindex >= symmetries_data[IPR].size()) return 0;

    for(int i=0;i<symmetries_data[IPR][Nindex].size();i++) 
      if(symmetries_data[IPR][Nindex][i] != " C1") 
	sum += symmetry_count_data[IPR][Nindex][i];
    return sum;
  } else {
    for(int i=0;i<symmetries_data[IPR][Nindex].size();i++) 
      if(symmetries_data[IPR][Nindex][i] == sym) 
	return symmetry_count_data[IPR][Nindex][i];
  }
  return 0;
}

vector<string> IsomerDB::symmetries(int N, bool IPR){ 
  int Nindex = (N-(IPR?60:20))/2;
  //  printf("symmetries(%d,%d) : %d = ",N,IPR,Nindex);
  if(Nindex >= symmetries_data[IPR].size()) return vector<string>();
  //  cout << symmetries_data[IPR][Nindex] << endl;

  return symmetries_data[IPR][Nindex]; 
}
