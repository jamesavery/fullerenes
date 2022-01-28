#pragma once


typedef uint16_t node_t; 

const node_t cubic_neighbours[60*3]  = {4,15,1,0,12,2,1,9,3,2,5,4,3,8,0,3,11,6,5,22,7,6,20,8,7,18,4,2,14,10,9,26,11,10,25,5,1,17,13,12,30,14,13,29,9,0,19,16,15,34,17,16,33,12,8,21,19,18,37,15,7,24,21,20,38,18,6,25,23,22,42,24,23,40,20,11,28,22,10,29,27,26,45,28,27,44,25,14,32,26,13,33,31,30,48,32,31,47,29,17,36,30,16,37,35,34,51,36,35,50,33,19,39,34,21,41,39,38,53,37,24,43,41,40,54,38,23,44,43,42,55,40,28,46,42,27,47,46,45,57,44,32,49,45,31,50,49,48,58,47,36,52,48,35,53,52,51,59,50,39,54,51,41,56,53,43,57,56,55,59,54,46,58,55,49,59,57,52,56,58};    
const node_t next_on_face[60*3] = {8,16,2,15,13,3,12,10,4,9,6,0,5,18,1,2,25,7,11,23,8,22,21,4,20,19,3,1,29,11,14,27,5,26,22,3,0,33,14,17,31,9,30,26,2,4,37,17,19,35,12,34,30,1,7,38,15,21,34,0,6,40,18,24,39,8,5,28,24,25,43,20,42,41,7,10,44,6,9,32,28,29,46,25,45,42,11,13,47,10,12,36,32,33,49,29,48,45,14,16,50,13,15,39,36,37,52,33,51,48,17,18,53,16,20,54,37,41,51,19,23,55,38,43,53,21,22,46,40,44,56,24,27,57,23,26,49,44,47,55,28,31,58,27,30,52,47,50,57,32,35,59,31,34,54,50,53,58,36,38,56,35,40,59,39,42,58,54,57,52,41,45,59,43,48,56,46,51,55,49};
const node_t prev_on_face[60*3] = {3,19,12,4,17,9,0,14,5,1,11,8,2,7,15,4,10,22,3,25,20,5,24,18,6,21,0,3,13,26,2,29,25,9,28,6,2,16,30,1,33,29,12,32,10,1,18,34,0,37,33,15,36,13,4,20,37,8,39,16,8,23,38,7,41,19,7,11,42,6,44,40,22,43,21,5,27,23,11,14,45,10,47,44,26,46,22,9,31,27,14,17,48,13,50,47,30,49,26,12,35,31,17,19,51,16,53,50,34,52,30,15,38,35,18,40,53,21,54,34,20,42,54,24,56,39,24,28,55,23,57,41,25,45,43,28,32,57,27,58,42,29,48,46,32,36,58,31,59,45,33,51,49,36,39,59,35,56,48,37,41,52,38,55,51,40,46,59,43,58,53,44,49,56,47,52,55,50,54,57};
const uint8_t face_right[60*3]   = {6,6,5,6,6,5,6,6,5,6,6,5,6,6,5,6,5,6,5,6,6,6,5,6,5,6,6,6,5,6,5,6,6,6,5,6,6,5,6,5,6,6,6,5,6,6,5,6,5,6,6,6,5,6,5,6,6,6,5,6,6,6,5,6,6,5,5,6,6,6,5,6,5,6,6,6,6,5,5,6,6,6,5,6,5,6,6,6,6,5,5,6,6,6,5,6,5,6,6,6,6,5,5,6,6,6,5,6,5,6,6,6,6,5,6,5,6,5,6,6,5,6,6,6,5,6,6,6,5,6,6,5,5,6,6,6,6,5,6,6,5,5,6,6,6,6,5,6,6,5,5,6,6,6,6,5,6,6,5,5,6,6,6,6,5,6,5,6,5,6,6,6,5,6,6,5,6,6,5,6};
const IsomerspaceForcefield::device_real_t X[60*3]= {1.429499999964262 ,-0.016451511315326996 ,3.2626215937469403 ,2.6090467778467734 ,0.1443134700710105 ,2.4208448897017227 ,2.6090467779088855 ,1.506336731552086 ,1.9005981835150294 ,1.4295000000186788 ,2.1873484191734245 ,2.420844740645028 ,0.7004999999712004 ,1.2462135273903634 ,3.2626215016248366 ,0.7290000000370741 ,3.0916899715509403 ,1.611976947633543 ,-0.7290000000370747 ,3.0916899715509403 ,1.6119769476335437 ,-1.429500000018679 ,2.1873484191734245 ,2.420844740645028 ,-0.7004999999712004 ,1.2462135273903634 ,3.262621501624837 ,3.041979587095356 ,1.7562905099906647 ,0.5918226904683742 ,2.3129795870195404 ,2.697425401748231 ,-0.24995407067044448 ,1.1795467778745952 ,3.351813175688808 ,0.24995362932763143 ,3.0419795870201627 ,-0.9145081899734363 ,1.611977239919614 ,3.492526364894045 ,-0.6543849858985488 ,0.24995392161370214 ,3.492526364894142 ,0.6543904525652654 ,-0.24995392161374932 ,0.7289999999930615 ,-1.2297531019674404 ,3.262621682267563 ,1.1795467778280952 ,-2.3316531593810463 ,2.4208450703444497 ,2.3129795869830883 ,-2.17717322875536 ,1.6119773320417181 ,-1.4294999999642615 ,-0.016451511315327346 ,3.2626215937469407 ,-0.7289999999930613 ,-1.2297531019674401 ,3.262621682267563 ,-2.609046777908885 ,1.5063367315520853 ,1.90059818351503 ,-2.6090467778467734 ,0.1443134700710103 ,2.420844889701723 ,-1.179546777874596 ,3.3518131756888074 ,0.24995362932763152 ,-2.3129795870195413 ,2.6974254017482306 ,-0.24995407067044398 ,-3.041979587095356 ,1.7562905099906643 ,0.5918226904683744 ,-3.795811684059282e-16 ,3.5125781571124883 ,-0.5918230748765796 ,2.3129795869832366 ,2.177178695422148 ,-1.611977332041806 ,1.1795467778281394 ,2.3316586260477674 ,-2.420845070344468 ,-6.5801150217778e-16 ,3.0126703137766833 ,-1.9005985132144685 ,3.0419795870203417 ,0.9145136566401717 ,-1.611977239919702 ,3.0419795870953092 ,-1.7562850433239707 ,-0.5918226904684047 ,2.6090467779089295 ,-1.5063312648854443 ,-1.900598183515047 ,2.6090467778468427 ,-0.1443080034043318 ,-2.420844889701742 ,2.3129795870194494 ,-2.6974199350815042 ,0.2499540706703975 ,7.639750375082463e-18 ,-3.0126648471099657 ,1.9005985132144512 ,4.2212170859923525e-16 ,-3.512572690445875 ,0.5918230748765495 ,1.1795467778746014 ,-3.35180770902225 ,-0.24995362932767815 ,-1.1795467778280946 ,-2.331653159381047 ,2.42084507034445 ,-3.0419795870201622 ,-0.9145081899734369 ,1.6119772399196148 ,-2.3129795869830883 ,-2.17717322875536 ,1.6119773320417183 ,-3.4925263648941423 ,0.6543904525652648 ,-0.24995392161374863 ,-3.492526364894045 ,-0.6543849858985491 ,0.24995392161370247 ,-2.3129795869832375 ,2.1771786954221475 ,-1.6119773320418056 ,-3.041979587020342 ,0.914513656640171 ,-1.6119772399197014 ,-1.1795467778281403 ,2.331658626047767 ,-2.420845070344468 ,0.728999999993055 ,1.2297585686340886 ,-3.262621682267386 ,-0.7289999999930559 ,1.2297585686340884 ,-3.262621682267386 ,1.4294999999642426 ,0.016456977981998015 ,-3.2626215937467635 ,1.4295000000187033 ,-2.1873429525068238 ,-2.4208447406450464 ,0.700499999971187 ,-1.2462080607236823 ,-3.26262150162466 ,0.7290000000371051 ,-3.091684504884463 ,-1.6119769476336312 ,-1.1795467778746012 ,-3.3518077090222493 ,-0.2499536293276779 ,-0.7290000000371049 ,-3.091684504884463 ,-1.6119769476336312 ,-2.312979587019449 ,-2.6974199350815047 ,0.24995407067039785 ,-3.0419795870953092 ,-1.7562850433239707 ,-0.5918226904684042 ,-2.609046777846843 ,-0.14430800340433222 ,-2.420844889701741 ,-2.60904677790893 ,-1.5063312648854446 ,-1.900598183515047 ,-1.4294999999642433 ,0.016456977981997793 ,-3.2626215937467635 ,-0.7004999999711874 ,-1.2462080607236825 ,-3.26262150162466 ,-1.4295000000187035 ,-2.1873429525068238 ,-2.4208447406450464};
