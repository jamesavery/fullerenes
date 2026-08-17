#include "fullerenes/config.hh"
#include "fullerenes/isomerdb.hh"

#include <stdexcept>
#include <cmath>
#include <cstring>

string IsomerDB::database_path = FULLERENE_DATABASE_PATH;

// A missing or corrupt database file is an environment/configuration error:
// there is no valid IsomerDB to return, and both prior behaviors -- abort()
// and returning an IsomerDB(-1) sentinel that unchecked callers dereference
// -- lose the diagnosis.  All readers fail loud through here instead.
[[noreturn]] static void throw_db_error(const string& where, const string& what,
                                        const string& filename)
{
  throw std::runtime_error(
      "IsomerDB::" + where + ": " + what + " '" + filename +
      "' (database root '" + IsomerDB::database_path +
      "'; set the FULLERENE_DATABASE_PATH CMake variable or assign "
      "IsomerDB::database_path)");
}

// Metadata for Cn isomers
vector<size_t> IsomerDB::Nisomers_data[2]                 = {
  {1,0,1,1,2,3,6,6,15,17,40,45,89,116,199,271,437,580,924,1205,1812,2385,3465,4478,6332,8149,11190,14246,19151,24109,31924,39718,51592,63761,81738,99918,126409,153493,191839,231017,285914,341658,419013,497529,604217,713319,860161,1008444,1207119,1408553,1674171,1942929,2295721,2650866,3114236,3580637,4182071,4787715,5566949,6344698,7341204,8339033,9604411,10867631,12469092,14059174,16066025,18060979,20558767,23037594,26142839,29202543,33022573,36798433,41478344,46088157,51809031,57417264,64353269,71163452,79538751,87738311,97841183,107679717,119761075,131561744,145976674,159999462,177175687,193814658,214127742,233846463, 257815889, 281006325, 309273526, 336500830, 369580714, 401535955, 440216206, 477420176, 522599564, 565900181, 618309598, 668662698, 729414880, 787556069, 857934016, 925042498, 1006016526, 1083451816, 1176632247, 1265323971, 1372440782, 1474111053, 1596482232, 1712934069, 1852762875, 1985250572, 2144943655, 2295793276, 2477017558, 2648697036, 2854536850, 3048609900, 3282202941, 3501931260, 3765465341, 4014007928, 4311652376, 4591045471, 4926987377, 5241548270, 5618445787, 5972426835, 6395981131, 6791769082, 7267283603, 7710782991, 8241719706, 8738236515, 9332065811, 9884604767, 10548218751, 11164542762, 11902015724, 12588998862, 13410330482, 14171344797, 15085164571, 15930619304, 16942010457, 17880232383, 19002055537, 20037346408, 21280571390, 22426253115, 23796620378, 25063227406, 26577912084, 27970034826, 29642262229, 31177474996, 33014225318, 34705254287, 36728266430, 38580626759, 40806395661, 42842199753, 45278616586, 47513679057, 50189039868, 52628839448, 55562506886, 58236270451, 61437700788, 64363670678, 67868149215, 71052718441, 74884539987, 78364039771, 82532990559, 86329680991, 90881152117, 95001297565, 99963147805, 104453597992, 109837310021, 114722988623, 120585261143, 125873325588, 132247999328},
  {1,0,0,0,0,1,1,1,2,5,7,9,24,19,35,46,86,134,187,259,450,616,823,1233,1799,2355,3342,4468,6063,8148,10774,13977,18769,23589,30683,39393,49878,62372,79362,98541,121354,151201,186611,225245,277930,335569,404667,489646,586264,697720,836497,989495,1170157,1382953,1628029,1902265,2234133,2601868,3024383,3516365,4071832,4690880,5424777,6229550,7144091,8187581,9364975,10659863,12163298,13809901,15655672,17749388, 20070486, 22606939, 25536557, 28700677, 32230861, 36173081, 40536922, 45278722, 50651799, 56463948, 62887775, 69995887, 77831323, 86238206, 95758929, 105965373, 117166528, 129476607, 142960479, 157402781, 173577766, 190809628, 209715141, 230272559, 252745513, 276599787, 303235792, 331516984, 362302637, 395600325, 431894257, 470256444, 512858451, 557745670, 606668511, 659140287, 716217922, 776165188, 842498881, 912274540, 987874095, 1068507788, 1156161307, 1247686189, 1348832364, 1454359806, 1568768524, 1690214836, 1821766896, 1958581588, 2109271290, 2266138871, 2435848971, 2614544391, 2808510141, 3009120113, 3229731630, 3458148016, 3704939275, 3964153268, 4244706701, 4533465777, 4850870260, 5178120469, 5531727283, 5900369830, 6299880577, 6709574675, 7158963073, 7620446934, 8118481242, 8636262789, 9196920285, 9768511147, 10396040696,11037658075,11730538496,12446446419,13221751502,14010515381,14874753568,15754940959,16705334454,17683643273,18744292915,19816289281,20992425825,22186413139,23475079272,24795898388,26227197453,27670862550,29254036711,30852950986,32581366295,34345173894,36259212641,38179777473,40286153024}
};

vector< vector<string> > IsomerDB::symmetries_data[2]     = {{{" Ih"},{},{"D6d"},{"D3h"},{" D2"," Td"},{"C2v","D5h"},{" C2"," D2"," D3","D3d","D3h"},{" C2"," Cs","C3v"},{" C1"," C2"," Cs"," D2","C2v","D2d","D3h","D6h"},{" C1"," C2"," D3","C2v","C3v","D3h"},{" C1"," C2"," C3"," Cs"," D2"," Td","C2v","C3v","D2h","D5d"},{" C1"," C2"," Cs"," D3","C2v"},{"  T"," C1"," C2"," Cs"," D2"," D3"," S4","C2v","D3d","D3h"},{" C1"," C2"," C3"," Cs","C2v"},{" C1"," C2"," Cs"," D2"," D3","C2h","C2v","D2h","D6d"},{" C1"," C2"," C3"," Cs"," D3","C2v","C3v","D3h","D5h"},{"  T"," C1"," C2"," C3"," Cs"," D2","C2h","C2v","C3v","D2d","D2h"},{" C1"," C2"," Cs"," D3","C2v","D3h"},{" C1"," C2"," C3"," Ci"," Cs"," D2"," D3"," Td","C2h","C2v","D2d","D2h","D3d"},{" C1"," C2"," C3"," Cs","C2v","C3v"},{" C1"," C2"," Cs"," D2"," D3"," D5"," Ih"," S4","C2h","C2v","C3v","D2d","D2h","D5d","D6h"},{" C1"," C2"," C3"," Cs"," D3","C2v","C3h","C3v","D3h"},{" C1"," C2"," C3"," Ci"," Cs"," D2","C2h","C2v","C3v","D2d","D2h"},{" C1"," C2"," Cs"," D3","C2v","C3v"},{"  T"," C1"," C2"," C3"," Cs"," D2"," D3"," S4"," S6"," Td","C2h","C2v","C3h","D2d","D2h","D3d","D3h"},{" C1"," C2"," C3"," Cs","C2v","C3v","D5h"},{" C1"," C2"," C3"," Ci"," Cs"," D2"," D3"," D6","C2h","C2v","C3v","D2d","D2h","D6d"},{" C1"," C2"," C3"," Cs"," D3","C2v","C3h","C3v","D3h"},{"  T"," C1"," C2"," C3"," Ci"," Cs"," D2"," S4"," Td","C2h","C2v","C3v","D2d"},{" C1"," C2"," C3"," Cs"," D3","C2v","D3h"},{" C1"," C2"," C3"," Ci"," Cs"," D2"," D3"," Ih"," S4","C2h","C2v","C3v","D3d","D3h","D5d","D5h"},{" C1"," C2"," C3"," Cs","C2v","C3v"},{" C1"," C2"," C3"," Ci"," Cs"," D2"," D3"," S4"," Td","C2h","C2v","C3v","D2d","D2h","D3d","D3h","D6h"},{" C1"," C2"," C3"," Cs"," D3","C2v","C3v","D3h"},{"  T"," C1"," C2"," C3"," Ci"," Cs"," D2"," S4","C2h","C2v","C3v","D2d","D2h"},{" C1"," C2"," C3"," Cs"," D3"," D5","C2v","D3h","D5h"},{"  T"," C1"," C2"," C3"," Ci"," Cs"," D2"," D3"," S4"," Td"," Th","C2h","C2v","C3h","C3v","D2d","D2h","D3d","D3h"},{" C1"," C2"," C3"," Cs","C2v","C3v"},{" C1"," C2"," C3"," Ci"," Cs"," D2"," D3","C2h","C2v","C3v","D2d","D2h","D3d","D3h","D6d","D6h"},{" C1"," C2"," C3"," Cs"," D3","C2v","C3h","C3v","D3h"},{"  T"," C1"," C2"," C3"," Ci"," Cs"," D2"," D5"," S4"," Td","C2h","C2v","C3v","D2d","D2h","D5d"},{" C1"," C2"," C3"," Cs"," D3","C2v","C3v"},{"  T"," C1"," C2"," C3"," Ci"," Cs"," D2"," D3"," S4","C2h","C2v","C3v","D2d","D2h","D3d","D3h"},{" C1"," C2"," C3"," Cs","C2v","C3v"},{" C1"," C2"," C3"," Ci"," Cs"," D2"," D3"," D6"," S4","C2h","C2v","C3h","C3v","D2d","D2h","D3d","D3h","D6h"},{" C1"," C2"," C3"," Cs"," D3"," D5","C2v","C3v","D3h","D5h"},{" C1"," C2"," C3"," Ci"," Cs"," D2"," S4"," Td","C2h","C2v","C3v","D2d","D2h"},{" C1"," C2"," C3"," Cs"," D3","C2v","C3v","D3h"},{"  T"," C1"," C2"," C3"," Ci"," Cs"," D2"," D3"," S4"," S6"," Th","C2h","C2v","C3h","C3v","D2d","D2h","D3d","D3h"},{" C1"," C2"," C3"," Cs","C2v","C3v"},{" C1"," C2"," C3"," Ci"," Cs"," D2"," D3"," D5"," D6"," S4"," Td","C2h","C2v","C3v","D2d","D2h","D5d","D5h","D6d"}},

 {{" Ih"},{},{},{},{},{"D5h"},{"D6d"},{"D3h"},{" D2"," Td"},{" D3","C2v","D3h"},{" D2"," D3"," Ih","C2v","D5d","D5h"},{" C2"," Cs","C2v","C3v"},{" C1"," C2"," Cs"," D2"," Td","C2v","D2d","D3d","D6h"},{" C1"," C2"," C3"," Cs"," D3","C2v"},{"  T"," C1"," C2"," Cs"," D2","C2v"},{" C1"," C2"," Cs","C2v","D5h"},{"  T"," C1"," C2"," C3"," Cs"," D2"," D3","C2v","D2h"},{" C1"," C2"," C3"," Cs","C2v","C3v"},{" C1"," C2"," Cs"," D2"," D3","C2v","C3v","D2d","D2h","D3d","D3h","D6d","D6h"},{" C1"," C2"," C3"," Cs"," D3","C2v"},{"  T"," C1"," C2"," C3"," Cs"," D2"," D5","C2v","D2d","D5d"},{" C1"," C2"," Cs"," D3","C2v","C3v"},{" C1"," C2"," C3"," Cs"," D2"," D3","C2v","D2d","D2h","D3d","D3h"},{" C1"," C2"," C3"," Cs","C2v","C3v"},{" C1"," C2"," Cs"," D2"," D3"," S4","C2h","C2v","C3v","D2d","D2h","D3d","D3h","D6h"},{" C1"," C2"," C3"," Cs"," D3"," D5","C2v","D5h"},{" C1"," C2"," C3"," Cs"," D2"," Td","C2h","C2v","C3v","D2h"},{" C1"," C2"," C3"," Cs"," D3","C2v","C3v","D3h"},{"  T"," C1"," C2"," C3"," Cs"," D2"," D3"," Th","C2h","C2v","C3h","D2d","D2h"},{" C1"," C2"," C3"," Cs","C2v","C3v"},{" C1"," C2"," C3"," Ci"," Cs"," D2"," D3"," D5"," D6"," Td","C2h","C2v","C3v","D2d","D2h","D5d","D5h","D6d"},{" C1"," C2"," C3"," Cs"," D3","C2v"},{"  T"," C1"," C2"," C3"," Ci"," Cs"," D2"," S4","C2h","C2v","C3v","D2d","D2h"},{" C1"," C2"," C3"," Cs"," D3","C2v"},{" C1"," C2"," C3"," Ci"," Cs"," D2"," D3"," S6","C2h","C2v","D2d","D3d","D3h"},{" C1"," C2"," C3"," Cs"," D5","C2v","C3v","D5h"},{"  T"," C1"," C2"," C3"," Cs"," D2"," D3"," D6"," S4","C2h","C2v","C3v","D2d","D2h","D3d","D3h","D6h"},{" C1"," C2"," C3"," Cs"," D3","C2v","C3h","C3v","D3h"},{"  T"," C1"," C2"," C3"," Ci"," Cs"," D2"," S4","C2h","C2v","C3v","D2d","D2h"},{" C1"," C2"," C3"," Cs"," D3","C2v","C3v"},{"  I","  T"," C1"," C2"," C3"," Ci"," Cs"," D2"," D3"," D5"," S4","C2h","C2v","C3h","D5d"},{" C1"," C2"," C3"," Cs","C2v","C3v"},{" C1"," C2"," C3"," Ci"," Cs"," D2"," D3"," D6"," S4","C2h","C2v","C3v","D2d","D2h","D3d","D3h","D6d","D6h"},{" C1"," C2"," C3"," Cs"," D3","C2v","D3h"},{"  T"," C1"," C2"," C3"," Ci"," Cs"," D2"," S4","C2h","C2v","C3v","D2d","D2h"},{" C1"," C2"," C3"," Cs"," D3"," D5","C2v","C3v","D3h","D5h"},{"  T"," C1"," C2"," C3"," Ci"," Cs"," D2"," D3"," Th","C2h","C2v","C3v","D2d","D3d","D3h"},{" C1"," C2"," C3"," Cs","C2v","C3v"},{"  T"," C1"," C2"," C3"," Ci"," Cs"," D2"," D3"," D6"," S4","C2h","C2v","C3v","D2d","D2h","D3d","D3h","D6h"},{" C1"," C2"," C3"," Cs"," D3","C2v","C3h"},{" C1"," C2"," C3"," Ci"," Cs"," D2"," D5"," S4"," Td","C2h","C2v","C3v","D2d","D2h","D5d","D5h"},{" C1"," C2"," C3"," Cs"," D3","C2v","C3v","D3h"},{"  T"," C1"," C2"," C3"," Ci"," Cs"," D2"," D3"," S4"," S6"," Td","C2h","C2v","C3h","D2d","D2h","D3d"},{" C1"," C2"," C3"," Cs","C2v","C3v"},{" C1"," C2"," C3"," Ci"," Cs"," D2"," D3"," D6"," S4"," Td","C2h","C2v","C3v","D2d","D2h","D3d","D6d"},{" C1"," C2"," C3"," Cs"," D3"," D5","C2v","C3v","D3h","D5h"},{"  T"," C1"," C2"," C3"," Ci"," Cs"," D2"," S4"," Td","C2h","C2v","D2d","D2h"},{" C1"," C2"," C3"," Cs"," D3","C2v","C3v","D3h"},{"  T"," C1"," C2"," C3"," Ci"," Cs"," D2"," D3"," S4","C2h","C2v","C3h","C3v","D2d","D2h","D3d","D3h"},{" C1"," C2"," C3"," Cs","C2v","C3v"},{" C1"," C2"," C3"," Ci"," Cs"," D2"," D3"," D5"," D6"," Ih"," S4","C2h","C2v","C3v","D2d","D2h","D3d","D3h","D5d","D6h"},{" C1"," C2"," C3"," Cs"," D3","C2v","C3h"},{" C1"," C2"," C3"," Ci"," Cs"," D2"," S4","C2h","C2v","C3v","D2d","D2h"},{" C1"," C2"," C3"," Cs"," D3","C2v","C3h","C3v","D3h"},{"  T"," C1"," C2"," C3"," Ci"," Cs"," D2"," D3"," S4"," S6"," Th","C2h","C2v","C3v","D2d","D2h","D3d","D3h"},{" C1"," C2"," C3"," Cs"," D5","C2v","C3v","D5h"},{" C1"," C2"," C3"," Ci"," Cs"," D2"," D3"," D6"," S4","C2h","C2v","C3v","D2d","D2h","D3d","D6d","D6h"},{" C1"," C2"," C3"," Cs"," D3","C2v","C3h","C3v","D3h"},{"  T"," C1"," C2"," C3"," Ci"," Cs"," D2"," S4","C2h","C2v","C3v","D2d","D2h"},{" C1"," C2"," C3"," Cs"," D3","C2v","C3v","D3h"},{"  T"," C1"," C2"," C3"," Ci"," Cs"," D2"," D3"," D5"," S4"," Th","C2h","C2v","C3v","D2d","D2h","D3d","D3h","D5d","D5h"}}
};

vector< vector<size_t> > IsomerDB::symmetry_count_data[2] = {{{1},{},{1},{1},{1,1},{2,1},{2,1,1,1,1},{3,2,1},{2,4,2,2,1,2,1,1},{7,5,1,2,1,1},{8,14,1,7,3,1,2,1,1,2},{23,11,6,1,4},{1,42,22,7,6,2,1,3,3,2},{69,22,2,19,4},{117,52,16,5,1,1,3,2,2},{195,37,2,25,2,6,1,1,2},{1,307,78,3,26,9,1,3,2,5,2},{470,62,38,1,8,1},{700,135,3,1,49,10,6,1,2,13,1,1,2},{1037,98,4,58,6,2},{1508,189,67,19,3,1,1,2,4,9,1,4,1,1,2},{2135,142,4,80,4,16,1,1,2},{2990,316,8,2,118,17,4,4,4,1,1},{4134,211,112,2,18,1},{1,5714,411,5,122,28,10,2,1,1,7,21,1,2,2,3,1},{7634,300,8,186,14,5,2},{10304,619,1,3,190,24,3,1,7,26,1,4,5,2},{13557,414,9,237,6,18,1,1,3},{1,18005,800,14,2,246,45,4,1,11,14,5,3},{23197,557,2,312,3,35,3},{30280,1146,15,5,371,39,12,1,1,16,28,2,2,2,2,2},{38548,742,15,380,28,5},{49590,1436,1,9,434,59,6,3,1,4,29,1,9,6,1,1,2},{62212,976,15,505,9,36,5,3},{1,79033,1945,24,10,596,52,1,16,43,8,4,5},{97936,1266,3,655,4,1,50,1,2},{1,123141,2412,20,12,646,80,19,5,1,1,13,38,2,2,4,4,5,3},{150939,1603,26,879,42,4},{187505,3200,4,20,972,70,9,16,28,2,3,3,1,1,3,2},{227934,2029,27,952,13,58,1,1,2},{2,280730,3801,40,14,1093,114,2,9,1,28,66,5,5,2,2},{337808,2542,3,1228,6,67,4},{1,412339,4954,37,29,1413,89,24,1,23,82,2,7,7,4,1},{492768,3126,38,1541,48,8},{596532,5872,5,26,1501,145,10,1,9,29,65,1,1,10,3,2,3,2},{707444,3845,42,1872,14,1,90,5,4,2},{850295,7403,59,41,2147,116,1,1,28,54,8,2,6},{1001569,4684,7,2079,9,88,5,3},{2,1195728,8713,49,44,2238,169,34,7,1,1,22,76,4,2,12,10,5,2},{1400184,5610,54,2609,84,12},{1660007,10787,10,52,2921,157,11,2,2,3,1,51,139,5,8,8,4,1,2}},

{{1},{},{},{},{},{1},{1},{1},{1,1},{1,2,2},{1,1,1,2,1,1},{3,3,1,2},{1,5,5,4,1,4,2,1,1},{6,6,1,3,1,2},{1,11,7,11,2,3},{16,16,6,7,1},{1,38,26,1,8,4,5,2,1},{89,26,3,13,2,1},{108,43,14,8,3,3,1,1,1,1,1,2,1},{169,49,3,30,3,5},{1,336,62,3,31,9,1,5,1,1},{488,73,44,2,8,1},{644,123,1,26,11,4,8,2,2,1,1},{1054,105,9,54,8,3},{1479,201,72,21,2,1,3,9,1,5,1,1,2,1},{2111,168,4,56,6,1,8,1},{2950,250,6,105,20,1,2,5,2,1},{4089,258,1,94,4,18,3,1},{1,5508,386,7,112,22,10,1,2,7,2,4,1},{7670,303,11,148,13,3},{9904,612,1,1,180,29,4,1,1,1,5,25,1,1,3,3,1,1},{13295,472,12,178,7,13},{2,17751,735,15,1,200,41,3,5,10,2,3,1},{22686,625,2,242,5,29},{29275,1052,12,3,269,36,11,1,4,13,1,2,4},{38285,773,24,274,1,29,6,1},{1,48013,1369,2,346,71,8,1,6,10,35,2,3,2,5,3,1},{60767,1117,13,433,13,26,1,1,1},{1,77072,1740,23,4,408,60,2,12,32,2,3,3},{96598,1382,6,517,5,31,2},{1,2,118321,2335,24,8,535,71,17,1,4,7,26,1,1},{148844,1676,31,612,35,3},{182804,2990,4,3,665,80,8,1,1,9,27,1,5,7,1,1,3,1},{222281,2171,33,710,12,37,1},{1,273339,3581,49,11,790,98,1,9,34,3,10,4},{331710,2792,7,982,10,2,57,4,2,3},{1,398804,4681,26,15,964,91,24,1,15,39,1,2,2,1},{485184,3275,44,1099,37,7},{1,579044,5752,8,17,1203,145,10,1,5,15,41,3,10,5,1,2,1},{692268,4051,48,1280,20,51,2},{828028,6732,59,19,1461,127,1,1,1,15,39,5,3,2,3,1},{982795,5083,7,1520,10,74,3,3},{1,1159784,8415,44,21,1561,176,34,11,2,1,26,70,2,3,4,2},{1375066,5852,68,1901,58,8},{1615320,10304,15,33,2039,163,19,2,2,1,20,89,5,6,6,4,1},{1893005,7049,62,2065,19,1,60,1,2,1},{3,2219523,11833,79,30,2374,207,9,1,29,41,3,1},{2590514,8685,15,2544,12,93,3,2},{1,3006861,14377,77,41,2704,182,36,1,20,61,1,1,6,7,3,4},{3503608,9817,104,2736,87,13},{4051036,17134,12,50,3115,286,15,5,2,1,12,37,96,4,12,5,2,4,1,3},{4675522,11800,73,3357,32,93,3},{5359303,18401,90,55,3637,179,3,27,90,6,4,5},{6211666,13865,26,3874,16,94,1,5,3},{2,7116261,23225,100,64,3946,302,50,11,1,1,35,80,1,7,3,1,1},{8167004,15775,119,4567,2,110,3,1},{9331690,27698,28,91,5028,300,20,1,2,32,61,6,6,7,1,3,1},{10636469,18521,118,4586,30,133,1,2,3},{4,12126580,30655,150,78,5232,383,16,47,138,7,6,2},{13781897,21759,22,6068,18,132,2,3},{1,15612518,36317,126,97,5982,330,51,3,8,1,58,163,2,3,4,2,1,3,2}}
};



// ------------------------------ Functions ------------------------------

// ---- Fortran formatted-field primitives (edit descriptors Aw, Iw, Fw.d) ----
//
// The database is written and read by Fortran FORMAT statements, so the
// text layer below is a column-exact emulation of the edit descriptors those
// statements use -- nothing more.  Input: Fortran BLANK='NULL' semantics
// (blanks anywhere in a numeric field are ignored, an all-blank field is 0).
// Output: right-justified in the field; a value that does not fit throws
// instead of the Fortran behaviour (a field of asterisks -- silent loss).

static string_view fortran_field(const string& s, int& pos, int len, const char* descriptor)
{
  if(pos < 0 || len < 0 || size_t(pos) + size_t(len) > s.size())
    throw std::invalid_argument(string("Fortran ") + descriptor + " field at column " +
                                to_string(pos+1) + " runs past the end of the record");
  string_view f(s.data() + pos, len);
  pos += len;
  return f;
}

void IsomerDB::fortran_readA(char *result, const string& s, int& pos, int len)
{
  string_view f = fortran_field(s, pos, len, "A");
  memcpy(result, f.data(), len);
}

long IsomerDB::fortran_readI(const string& s, int& pos, int len)
{
  string_view f = fortran_field(s, pos, len, "I");
  long v = 0; bool negative = false, sign_seen = false, digit_seen = false;
  for(char c: f){
    if(c == ' ') continue;
    if((c == '-' || c == '+') && !sign_seen && !digit_seen){ sign_seen = true; negative = (c=='-'); continue; }
    if(c >= '0' && c <= '9'){ v = 10*v + (c-'0'); digit_seen = true; continue; }
    throw std::invalid_argument("Fortran I" + to_string(len) + " field '" + string(f) + "' is not an integer");
  }
  return negative? -v : v;
}

double IsomerDB::fortran_readF(const string& s, int& pos, int len)
{
  string_view f = fortran_field(s, pos, len, "F");
  string compact;                          // the field with blanks removed
  for(char c: f) if(c != ' ') compact += c;
  if(compact.empty()) return 0;
  // Fortran F input accepts [sign] digits [. digits] [exponent]; the database
  // only ever holds plain decimals, and a field without a decimal point would
  // silently take an implied scale (Fw.d) -- refuse anything but sign, digits
  // and one point.
  size_t points = 0, digits = 0;
  for(size_t i=0;i<compact.size();i++){
    char c = compact[i];
    if(c == '.') points++;
    else if(c >= '0' && c <= '9') digits++;
    else if(!((c == '-' || c == '+') && i == 0))
      throw std::invalid_argument("Fortran F" + to_string(len) + " field '" + string(f) + "' is not a decimal number");
  }
  if(points != 1 || digits == 0)
    throw std::invalid_argument("Fortran F" + to_string(len) + " field '" + string(f) + "' is not a decimal number");
  return strtod(compact.c_str(), nullptr);
}

string IsomerDB::fortran_writeI(long v, int len)
{
  string digits = to_string(v);
  if((int)digits.size() > len)
    throw std::invalid_argument("value " + digits + " does not fit a Fortran I" + to_string(len) + " field");
  return string(len - digits.size(), ' ') + digits;
}

string IsomerDB::fortran_writeF(double v, int len, int decimals)
{
  char buf[64];
  int n = snprintf(buf, sizeof buf, "%*.*f", len, decimals, v);
  if(n < 0 || n > len)
    throw std::invalid_argument("value " + to_string(v) + " does not fit a Fortran F" +
                                to_string(len) + "." + to_string(decimals) + " field");
  return string(buf, n);
}

// ---- One text record: cursor with the file/line context for diagnostics ----
namespace {

struct PDBRecord {
  const string& line; const string& filename; size_t lineno; int pos = 0;

  [[noreturn]] void fail(const string& what) const {
    throw_db_error("readPDB", "line " + to_string(lineno) + ": " + what +
                   " in record '" + line + "' of database file", filename);
  }
  string A(int w)     { try { string b(w, ' '); IsomerDB::fortran_readA(b.data(), line, pos, w); return b; } catch(const std::invalid_argument& e){ fail(e.what()); } }
  long   I(int w)     { try { return IsomerDB::fortran_readI(line, pos, w); } catch(const std::invalid_argument& e){ fail(e.what()); } }
  double F(int w)     { try { return IsomerDB::fortran_readF(line, pos, w); } catch(const std::invalid_argument& e){ fail(e.what()); } }
  // wX skip, checking that the skipped columns hold what the writing FORMAT
  // put there (a misaligned parse cannot pass this).
  void   X(int w, string_view expect) {
    string_view got;
    try { got = fortran_field(line, pos, w, "X"); } catch(const std::invalid_argument& e){ fail(e.what()); }
    if(got != expect) fail("expected '" + string(expect) + "' at column " + to_string(pos-w+1) + ", found '" + string(got) + "'");
  }
  void   done() const { if(size_t(pos) != line.size()) fail("trailing characters after column " + to_string(pos)); }
  // Values stored in one-byte fields must fit them (silent truncation is a bug).
  u_int8_t u8(long v, const char* what) const {
    if(v < 0 || v > 255) fail(string(what) + " = " + to_string(v) + " does not fit the 8-bit database field");
    return u_int8_t(v);
  }
};

// Structural checks shared by both text layouts.
void check_rspi(const PDBRecord& r, const u_int8_t* RSPI, int N)
{
  const int Nfaces = N/2 + 2;
  for(int i=0;i<12;i++){
    if(RSPI[i] < 1 || RSPI[i] > Nfaces)
      r.fail("RSPI[" + to_string(i) + "] = " + to_string(RSPI[i]) + " outside [1," + to_string(Nfaces) + "]");
    if(i > 0 && RSPI[i] <= RSPI[i-1])
      r.fail("RSPI not strictly ascending at position " + to_string(i));
  }
}

// Compact record: A3,12I3,{5I2,6I2 | 3I2},I2,I1,F7.5[,I7],6I3
IsomerDB::Entry parse_compact_record(PDBRecord& r, int N, bool IPR, bool with_ncycham)
{
  IsomerDB::Entry e{};
  string group = r.A(3); memcpy(e.group, group.data(), 3);
  for(int i=0;i<12;i++) e.RSPI[i] = r.u8(r.I(3), "RSPI");
  if(!IPR){
    for(int i=0;i<5;i++) e.PNI[i] = r.u8(r.I(2), "PNI");
    for(int i=0;i<6;i++) e.HNI[i] = r.u8(r.I(2), "HNI");
  } else {
    // Only HNI k=3..5 are stored: every pentagon of an IPR isomer has 0
    // pentagon neighbours and no hexagon has fewer than 3 hexagon neighbours.
    e.PNI[0] = 12;
    for(int i=3;i<6;i++) e.HNI[i] = r.u8(r.I(2), "HNI");
  }
  e.NeHOMO    = r.u8(r.I(2), "NeHOMO");
  e.NedgeHOMO = r.u8(r.I(1), "NedgeHOMO");
  e.HLgap     = r.F(7);
  if(with_ncycham){
    long h = r.I(7);
    if(h < 0 || h > INT32_MAX) r.fail("ncycham = " + to_string(h) + " out of range");
    e.ncycham = int(h);
  }
  for(int i=0;i<6;i++) e.INMR[i] = r.u8(r.I(3), "INMR");
  r.done();
  check_rspi(r, e.RSPI, N);
  return e;
}

// Verbose record (the program's printed isomer list), spiral.f/isomer.f
// format 608 (with Hamilton-cycle count) / 607 (without):
//   1X,I8,2X,A3,1X,12I4,2X,'(',5(I2,','),I2,')  ',I2,2X,'(',6(I2,','),I3,')  ',
//   F8.5,2X,I2,1X,I2,1X,F8.5,1X,A6[,2X,I9],2X,3(I3,' x',I3,:,',')
// i.e. number, group, RSPI, (PNI k=0..5), Np, (HNI k=0..6), sigma_h, NeHOMO,
// NedgeHOMO, HLgap, closed/open[, ncycham], NMR pattern.  The columns the
// compact format does not store are all derived (util.f IPentInd / HexInd,
// isomer.f Printdatabase) and are verified against the stored ones here.
// The leading isomer number is the record's index in the FULL list of the
// size; a filtered listing (the Nontrivial/ files) keeps the original
// numbers, so the invariant the reader can hold every file to is "strictly
// increasing" -- a full list's completeness is the caller's check against
// number_isomers().
IsomerDB::Entry parse_verbose_record(PDBRecord& r, int N, bool IPR, bool with_ncycham, long& last_number)
{
  IsomerDB::Entry e{};
  r.X(1, " ");
  long number = r.I(8);
  if(number <= last_number)
    r.fail("isomer number " + to_string(number) + " does not increase (previous " + to_string(last_number) + ")");
  last_number = number;
  r.X(2, "  ");
  string group = r.A(3); memcpy(e.group, group.data(), 3);
  r.X(1, " ");
  for(int i=0;i<12;i++) e.RSPI[i] = r.u8(r.I(4), "RSPI");

  long PNI[6], HNI[7];
  r.X(3, "  (");
  for(int k=0;k<6;k++){ PNI[k] = r.I(2); r.X(1, k<5? ",": ")"); }
  r.X(2, "  ");
  long Np = r.I(2);
  r.X(3, "  (");
  for(int k=0;k<6;k++){ HNI[k] = r.I(2); r.X(1, ","); }
  HNI[6] = r.I(3);
  r.X(1, ")");
  double sigmah = r.F(10);            // 2X,F8.5
  r.X(1, " ");
  long NeHOMO = r.I(3), NedgeHOMO = r.I(3);
  double HLgap = r.F(9);              // 1X,F8.5
  r.X(1, " ");
  string shell = r.A(6);
  long ncycham = 0;
  if(with_ncycham){
    r.X(1, " ");
    ncycham = r.I(10);                // 2X,I9 written; CompressDatabase reads 1X,I10
    if(ncycham <= 0) r.fail("Hamilton-cycle count " + to_string(ncycham) + " in a file declaring Hamilton counts");
    if(ncycham > INT32_MAX) r.fail("ncycham = " + to_string(ncycham) + " out of range");
  }
  r.X(2, "  ");
  // NMR pattern: 1..3 groups of I3,' x',I3 separated by ','.
  int nmr_len = int(r.line.size()) - r.pos;
  if(nmr_len != 8 && nmr_len != 17 && nmr_len != 26)
    r.fail("NMR pattern '" + r.line.substr(r.pos) + "' is not 1-3 'nnn xnnn' groups");
  for(int k=0;k<(nmr_len+1)/9;k++){
    if(k > 0) r.X(1, ",");
    e.INMR[2*k]   = r.u8(r.I(3), "INMR");
    r.X(2, " x");
    e.INMR[2*k+1] = r.u8(r.I(3), "INMR");
  }
  r.done();

  // The derived columns must agree with the stored ones.
  long pni_sum = 0, pni_moment = 0;
  for(int k=0;k<5;k++){ pni_sum += PNI[k]; pni_moment += k*PNI[k]; }
  pni_moment += 5*PNI[5];
  if(pni_sum + PNI[5] != 12) r.fail("pentagon neighbour indices do not sum to 12");
  if(Np != pni_moment/2)     r.fail("Np = " + to_string(Np) + " != IPentInd(PNI) = " + to_string(pni_moment/2));
  long hni_sum = 0;
  for(int k=0;k<7;k++) hni_sum += HNI[k];
  if(hni_sum != N/2 - 10)    r.fail("hexagon neighbour indices do not sum to " + to_string(N/2-10));
  // sigma_h = standard deviation of the hexagon-neighbour count over the
  // hexagons (util.f HexInd).  The C110-C120 listings were printed by a
  // program version whose HexInd summed only k=3..6, so both readings are
  // accepted -- the column is derived and dropped either way.
  auto hexind = [&](int kmin){
    long n = 0, hk = 0, hk2 = 0;
    for(int k=kmin;k<7;k++){ n += HNI[k]; hk += k*HNI[k]; hk2 += k*k*HNI[k]; }
    return n? sqrt(fabs(double(hk2)/n - pow(double(hk)/n,2))) : 0.0;
  };
  const double sigmah_all = hexind(0), sigmah_k3 = hexind(3), sigmah_tol = 5.1e-6;   // F8.5 rounding
  if(fabs(sigmah - sigmah_all) > sigmah_tol && fabs(sigmah - sigmah_k3) > sigmah_tol)
    r.fail("sigma_h = " + to_string(sigmah) + " != HexInd(HNI) = " + to_string(sigmah_all) +
           " (nor the legacy k>=3 value " + to_string(sigmah_k3) + ")");
  string shell_expected = (2*NedgeHOMO == NeHOMO)? "closed" : "open  ";
  if(shell != shell_expected) r.fail("shell '" + shell + "' inconsistent with NeHOMO/NedgeHOMO");
  if(IPR && (PNI[0] != 12 || HNI[0] || HNI[1] || HNI[2]))
    r.fail("non-IPR neighbour indices in an IPR file");

  for(int k=0;k<5;k++) e.PNI[k] = r.u8(PNI[k], "PNI");
  for(int k=0;k<6;k++) e.HNI[k] = r.u8(HNI[k], "HNI");
  e.NeHOMO    = r.u8(NeHOMO, "NeHOMO");
  e.NedgeHOMO = r.u8(NedgeHOMO, "NedgeHOMO");
  e.HLgap     = HLgap;
  e.ncycham   = int(ncycham);
  check_rspi(r, e.RSPI, N);
  return e;
}

// A verbose data record starts with the isomer number right-justified in
// columns 1-9 (1X,I8); the trailer lines the program prints afterwards
// ("Lowest number of Hamiltonian cycles ...", " Highest Np= ...") do not.
bool looks_like_verbose_record(const string& line)
{
  if(line.size() < 9) return false;
  bool digit = false;
  for(int i=0;i<9;i++){
    char c = line[i];
    if(c >= '0' && c <= '9') digit = true;
    else if(c != ' ') return false;
  }
  return digit;
}

void rtrim(string& s)
{
  while(!s.empty() && (s.back() == '\r' || s.back() == ' ' || s.back() == '\t')) s.pop_back();
}

}  // namespace

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
  fwrite(entries.data(),sizeof(Entry),Nisomers,f);

  bool ok = !ferror(f);
  fclose(f);
  return ok;
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
    return f.good();
  }

IsomerDB IsomerDB::readBinary(const string filename){
  FILE *f = fopen(filename.c_str(),"rb");
  if(!f)
    throw_db_error("readBinary", string("cannot open database file (")
                   + strerror(errno) + ")", filename);

  u_int16_t header;
  u_int32_t n;
  if(fread(&header,2,1,f) != 1 || fread(&n,4,1,f) != 1){
    fclose(f);
    throw_db_error("readBinary", "truncated header in database file", filename);
  }

  // Validate the declared entry count against the actual file size BEFORE
  // sizing any buffer by it: a corrupt count would otherwise force a huge
  // allocation or a silent short read.
  long data_start = ftell(f);
  fseek(f, 0, SEEK_END);
  long data_bytes = ftell(f) - data_start;
  fseek(f, data_start, SEEK_SET);
  if(data_bytes < 0 || (long long)n * (long long)sizeof(Entry) != (long long)data_bytes){
    fclose(f);
    throw_db_error("readBinary",
                   "corrupt database file (declares " + to_string(n) +
                   " entries of " + to_string(sizeof(Entry)) +
                   " bytes but holds " + to_string(data_bytes) +
                   " data bytes):", filename);
  }

  IsomerDB DB(header & 0xff, header >> 8 & 1, header >> 9 & 1);
  DB.Nisomers = (int)n;
  DB.entries.resize(n);
  if(n > 0 && fread(DB.entries.data(),sizeof(Entry),n,f) != n){
    fclose(f);
    throw_db_error("readBinary", "read error in database file", filename);
  }
  fclose(f);
  return DB;
}

IsomerDB::Entry IsomerDB::getIsomer(int N, int isomer, bool IPR){
  string filename = database_path+"/binary/c"+pad_string(to_string(N),3)+(IPR?"IPR":"all")+".bin";
  FILE *f = fopen(filename.c_str(),"rb");
  if(!f)
    throw_db_error("getIsomer", string("cannot open database file (")
                   + strerror(errno) + ")", filename);

  u_int16_t header;
  u_int32_t Nisomers;
  if(fread(&header,2,1,f) != 1 || fread(&Nisomers,4,1,f) != 1){
    fclose(f);
    throw_db_error("getIsomer", "truncated header in database file", filename);
  }
  if(isomer < 1 || (u_int32_t)isomer > Nisomers){
    fclose(f);
    throw std::out_of_range(
        "IsomerDB::getIsomer: isomer index " + to_string(isomer) +
        " outside [1, " + to_string(Nisomers) + "] for C" + to_string(N) +
        (IPR ? " (IPR)" : "") + " in '" + filename + "'");
  }

  Entry e;
  fseek(f,(long)(isomer-1)*sizeof(Entry),SEEK_CUR);
  if(fread(&e,sizeof(Entry),1,f) != 1){
    fclose(f);
    throw_db_error("getIsomer", "read error at isomer " + to_string(isomer)
                   + " in database file", filename);
  }
  fclose(f);
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

string IsomerDB::PDBfilename(int N, bool IPR, string extension)
{
  return database_path + (IPR? "/IPR/c" : "/All/c") + pad_string(to_string(N),3,'0')
                       + (IPR? "IPR" : "all") + extension + ".database";
}

IsomerDB IsomerDB::readPDB(int N, bool IPR, string extension)
{
  const string filename = PDBfilename(N, IPR, extension);
  IsomerDB DB = readPDB(filename);
  if(DB.N != N || DB.IPR != IPR)
    throw_db_error("readPDB", "header says C" + to_string(DB.N) + (DB.IPR? " IPR" : " all") +
                   " but C" + to_string(N) + (IPR? " IPR" : " all") + " was requested from database file",
                   filename);
  return DB;
}

IsomerDB IsomerDB::readPDB(const string& filename)
{
  ifstream dbfile(filename.c_str());
  if(!dbfile)
    throw_db_error("readPDB", string("cannot open database file (")
                   + strerror(errno) + ")", filename);

  // Header.  Compact files: Format(I3,2I1) -- exactly five characters.
  // Verbose listings: Format(I5,2I2) as printed by spiral.f (some hand-fixed
  // to "N IP IH"), read list-directed by util.f CompressDatabase -- three
  // blank-separated integers.
  string line;
  if(!getline(dbfile,line))
    throw_db_error("readPDB", "empty database file", filename);
  rtrim(line);
  vector<string> tokens;
  for(size_t i=0;i<line.size();){
    if(line[i] == ' '){ i++; continue; }
    size_t j = i; while(j < line.size() && line[j] != ' ') j++;
    tokens.push_back(line.substr(i,j-i)); i = j;
  }
  long Nread = 0, IP = -1, IH = -1;
  bool verbose = false;
  try {
    if(tokens.size() == 1 && line.size() == 5){
      int pos = 0;
      Nread = fortran_readI(line,pos,3); IP = fortran_readI(line,pos,1); IH = fortran_readI(line,pos,1);
    } else if(tokens.size() == 3){
      int pos = 0;
      Nread = fortran_readI(tokens[0],pos,tokens[0].size()); pos = 0;
      IP    = fortran_readI(tokens[1],pos,tokens[1].size()); pos = 0;
      IH    = fortran_readI(tokens[2],pos,tokens[2].size());
      verbose = true;
    }
  } catch(const std::invalid_argument&){ /* reported below */ }
  if(Nread < 20 || Nread == 22 || (Nread & 1) || (IP != 0 && IP != 1) || (IH != 0 && IH != 1))
    throw_db_error("readPDB", "malformed header line '" + line + "' in database file", filename);

  IsomerDB DB(int(Nread), IP == 1, IH == 1);
  size_t lineno = 1;
  bool in_trailer = false;
  long last_number = 0;
  while(getline(dbfile,line)){
    lineno++;
    rtrim(line);
    PDBRecord r{line, filename, lineno};
    if(!verbose){
      if(line.empty()) r.fail("blank line");
      DB.entries.push_back(parse_compact_record(r, DB.N, DB.IPR, DB.with_ncycham));
    } else {
      // Records, then the program's summary trailer; a record after the
      // trailer has begun means the file is scrambled, not merely trailed.
      if(!looks_like_verbose_record(line)){ in_trailer = true; continue; }
      if(in_trailer) r.fail("isomer record after the trailer lines");
      DB.entries.push_back(parse_verbose_record(r, DB.N, DB.IPR, DB.with_ncycham, last_number));
    }
  }

  DB.Nisomers = DB.entries.size();
  for(int i=0;i<DB.Nisomers;i++){
    const Entry &e(DB.entries[i]);
    DB.RSPIindex[vector<int>(e.RSPI,e.RSPI+12)] = i;
  }
  return DB;
}

// Compact text format, formats 1003 (header) and 1004/1007/1008/1009
// (records) of isomer.f -- the same edit descriptors readPDB parses, so a
// compact file round-trips byte for byte.
bool IsomerDB::writePDB(const string& filename) const
{
  ofstream f(filename.c_str());
  if(!f){
    cerr << "Couldn't open database file " << filename << " for writing: " << strerror(errno) << ".\n";
    return false;
  }
  f << fortran_writeI(N,3) << fortran_writeI(IPR,1) << fortran_writeI(with_ncycham,1) << "\n";
  string record;
  for(const Entry& e: entries){
    record.assign(e.group, 3);
    for(int i=0;i<12;i++) record += fortran_writeI(e.RSPI[i],3);
    if(!IPR){
      for(int i=0;i<5;i++) record += fortran_writeI(e.PNI[i],2);
      for(int i=0;i<6;i++) record += fortran_writeI(e.HNI[i],2);
    } else
      for(int i=3;i<6;i++) record += fortran_writeI(e.HNI[i],2);
    record += fortran_writeI(e.NeHOMO,2);
    record += fortran_writeI(e.NedgeHOMO,1);
    record += fortran_writeF(e.HLgap,7,5);
    if(with_ncycham) record += fortran_writeI(e.ncycham,7);
    for(int i=0;i<6;i++) record += fortran_writeI(e.INMR[i],3);
    f << record << "\n";
  }
  return f.good();
}

int64_t IsomerDB::number_isomers(int N, const string& sym, bool IPR){
  if(N < (IPR?60:20) || (N & 1)) return 0;   // below the family's base size, or odd
  int Nindex = (N-(IPR?60:20))/2;
  if(Nindex >= (int)Nisomers_data[IPR].size()) return 0;

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
  if(N < (IPR?60:20) || (N & 1)) return vector<string>();
  int Nindex = (N-(IPR?60:20))/2;
  if(Nindex >= (int)symmetries_data[IPR].size()) return vector<string>();

  return symmetries_data[IPR][Nindex];
}
