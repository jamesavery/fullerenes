#include <sys/stat.h>
#include <numeric>
#include <random>
#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/gpu/batch_queue.hh"
#include "fullerenes/gpu/isomer_batch.hh"
#include "fullerenes/gpu/kernels.hh"
#include "fullerenes/gpu/cuda_io.hh"
#include "fullerenes/gpu/cu_array.hh"
#include "fullerenes/progress_bar.hh"
const std::unordered_map<size_t,size_t> num_fullerenes = {{20,1},{22,0},{24,1},{26,1},{28,2},{30,3},{32,6},{34,6},{36,15},{38,17},{40,40},{42,45},{44,89},{46,116},{48,199},{50,271},{52,437},{54,580},{56,924},{58,1205},{60,1812},{62,2385},{64,3465},{66,4478},{68,6332},{70,8149},{72,11190},{74,14246},{76,19151},{78,24109},{80,31924},{82,39718},{84,51592},{86,63761},{88,81738},{90,99918},{92,126409},{94,153493},{96,191839},{98,231017},{100,285914},{102,341658},{104,419013},{106,497529},{108,604217},{110,713319},{112,860161},{114,1008444},{116,1207119},{118,1408553},{120,1674171},{122,1942929},{124,2295721},{126,2650866},{128,3114236},{130,3580637},{132,4182071},{134,4787715},{136,5566949},{138,6344698},{140,7341204},{142,8339033},{144,9604411},{146,10867631},{148,12469092},{150,14059174},{152,16066025},{154,18060979},{156,20558767},{158,23037594},{160,26142839},{162,29202543},{164,33022573},{166,36798433},{168,41478344},{170,46088157},{172,51809031},{174,57417264},{176,64353269},{178,71163452},{180,79538751},{182,87738311},{184,97841183},{186,107679717},{188,119761075},{190,131561744},{192,145976674},{194,159999462},{196,177175687},{198,193814658},{200,214127742},{202,233846463},{204,257815889},{206,281006325},{208,309273526},{210,336500830},{212,369580714},{214,401535955},{216,440216206},{218,477420176},{220,522599564},{222,565900181},{224,618309598},{226,668662698},{228,729414880},{230,787556069},{232,857934016},{234,925042498},{236,1006016526},{238,1083451816},{240,1176632247},{242,1265323971},{244,1372440782},{246,1474111053},{248,1596482232},{250,1712934069},{252,1852762875},{254,1985250572},{256,2144943655},{258,2295793276},{260,2477017558},{262,2648697036},{264,2854536850},{266,3048609900},{268,3282202941},{270,3501931260},{272,3765465341},{274,4014007928},{276,4311652376},{278,4591045471},{280,4926987377},{282,5241548270},{284,5618445787},{286,5972426835},{288,6395981131},{290,6791769082},{292,7267283603},{294,7710782991},{296,8241719706},{298,8738236515},{300,9332065811},{302,9884604767},{304,10548218751},{306,11164542762},{308,11902015724},{310,12588998862},{312,13410330482},{314,14171344797},{316,15085164571},{318,15930619304},{320,16942010457},{322,17880232383},{324,19002055537},{326,20037346408},{328,21280571390},{330,22426253115},{332,23796620378},{334,25063227406},{336,26577912084},{338,27970034826},{340,29642262229},{342,31177474996},{344,33014225318},{346,34705254287},{348,36728266430},{350,38580626759},{352,40806395661},{354,42842199753},{356,45278616586},{358,47513679057},{360,50189039868},{362,52628839448},{364,55562506886},{366,58236270451},{368,61437700788},{370,64363670678},{372,67868149215},{374,71052718441},{376,74884539987},{378,78364039771},{380,82532990559},{382,86329680991},{384,90881152117},{386,95001297565},{388,99963147805},{390,104453597992},{392,109837310021},{394,114722988623},{396,120585261143},{398,125873325588},{400,132247999328}};
const std::unordered_map<size_t,size_t> num_IPR_fullerenes = {{20,0},{22,0},{24,0},{26,0},{28,0},{30,0},{32,0},{34,0},{36,0},{38,0},{40,0},{42,0},{44,0},{46,0},{48,0},{50,0},{52,0},{54,0},{56,0},{58,0},{60,1},{62,0},{64,0},{66,0},{68,0},{70,1},{72,1},{74,1},{76,2},{78,5},{80,7},{82,9},{84,24},{86,19},{88,35},{90,46},{92,86},{94,134},{96,187},{98,259},{100,450},{102,616},{104,823},{106,1233},{108,1799},{110,2355},{112,3342},{114,4468},{116,6063},{118,8148},{120,10774},{122,13977},{124,18769},{126,23589},{128,30683},{130,39393},{132,49878},{134,62372},{136,79362},{138,98541},{140,121354},{142,151201},{144,186611},{146,225245},{148,277930},{150,335569},{152,404667},{154,489646},{156,586264},{158,697720},{160,836497},{162,989495},{164,1170157},{166,1382953},{168,1628029},{170,1902265},{172,2234133},{174,2601868},{176,3024383},{178,3516365},{180,4071832},{182,4690880},{184,5424777},{186,6229550},{188,7144091},{190,8187581},{192,9364975},{194,10659863},{196,12163298},{198,13809901},{200,15655672},{202,17749388},{204,20070486},{206,22606939},{208,25536557},{210,28700677},{212,32230861},{214,36173081},{216,40536922},{218,45278722},{220,50651799},{222,56463948},{224,62887775},{226,69995887},{228,77831323},{230,86238206},{232,95758929},{234,105965373},{236,117166528},{238,129476607},{240,142960479},{242,157402781},{244,173577766},{246,190809628},{248,209715141},{250,230272559},{252,252745513},{254,276599787},{256,303235792},{258,331516984},{260,362302637},{262,395600325},{264,431894257},{266,470256444},{268,512858451},{270,557745670},{272,606668511},{274,659140287},{276,716217922},{278,776165188},{280,842498881},{282,912274540},{284,987874095},{286,1068507788},{288,1156161307},{290,1247686189},{292,1348832364},{294,1454359806},{296,1568768524},{298,1690214836},{300,1821766896},{302,1958581588},{304,2109271290},{306,2266138871},{308,2435848971},{310,2614544391},{312,2808510141},{314,3009120113},{316,3229731630},{318,3458148016},{320,3704939275},{322,3964153268},{324,4244706701},{326,4533465777},{328,4850870260},{330,5178120469},{332,5531727283},{334,5900369830},{336,6299880577},{338,6709574675},{340,7158963073},{342,7620446934},{344,8118481242},{346,8636262789},{348,9196920285},{350,9768511147},{352,10396040696},{354,11037658075},{356,11730538496},{358,12446446419},{360,13221751502},{362,14010515381},{364,14874753568},{366,15754940959},{368,16705334454},{370,17683643273},{372,18744292915},{374,19816289281},{376,20992425825},{378,22186413139},{380,23475079272},{382,24795898388},{384,26227197453},{386,27670862550},{388,29254036711},{390,30852950986},{392,32581366295},{394,34345173894},{396,36259212641},{398,38179777473},{400,40286153024}};

int main(int ac, char **argv){
    
    int N_begin                = strtol(argv[1],0,0);     // Argument 1: Number of vertices N
    int N_limit                = strtol(argv[2],0,0);     // Argument 2: Number of vertices N
    int sample_size            = strtol(argv[3],0,0);     // Argument 3: smaple size
    std::string suffix         = argv[4];     // Argument 3: smaple size
    bool IPR = false;
    string path = "isomerspace_stats/" +  suffix + "/";
    mkdir("isomerspace_stats",0777);
    mkdir(path.data(), 0777);
    for(int N = N_begin; N < N_limit+1; N+=2){
        auto n_isomers = IPR ? num_IPR_fullerenes.find(N)->second : num_fullerenes.find(N)->second;
        auto batch_size = min(sample_size, (int)n_isomers);
        using namespace gpu_kernels;
        IsomerBatch d_validation(N, batch_size, DEVICE_BUFFER);
        IsomerBatch h_validation(N, batch_size, HOST_BUFFER);
        
        std::vector<int> random_IDs(n_isomers);
        std::iota(random_IDs.begin(), random_IDs.end(), 0);
        std::shuffle(random_IDs.begin(), random_IDs.end(), std::mt19937{42});
        
        std::vector<int> sorted_IDs(random_IDs.begin(), random_IDs.begin()+batch_size);
        std::sort(sorted_IDs.begin(), sorted_IDs.end());
        std::queue<int> ID_queue;
        for(int element : sorted_IDs){
            ID_queue.push(element);
        }

        BuckyGen::buckygen_queue Q = BuckyGen::start(N,IPR, false);
        FullereneDual FD;
        bool more_to_generate = true;
        

        cuda_io::IsomerQueue BQ(N);
        cuda_io::IsomerQueue BQ2(N);
        std::queue<std::tuple<Polyhedron, size_t, IsomerStatus>> out_queue;
        ProgressBar progress_bar = ProgressBar('#',30);

        int I = 0;
        while(more_to_generate)
        {
            more_to_generate &= BuckyGen::next_fullerene(Q,FD);
            if(!more_to_generate)break;
            if(I == ID_queue.front()){
                FD.update();
                PlanarGraph G = FD.dual_graph();
                //G.neighbours = {{4,15,1},{0,12,2},{1,9,3},{2,5,4},{3,8,0},{3,11,6},{5,22,7},{6,20,8},{7,18,4},{2,14,10},{9,26,11},{10,25,5},{1,17,13},{12,30,14},{13,29,9},{0,19,16},{15,34,17},{16,33,12},{8,21,19},{18,37,15},{7,24,21},{20,38,18},{6,25,23},{22,42,24},{23,40,20},{11,28,22},{10,29,27},{26,45,28},{27,44,25},{14,32,26},{13,33,31},{30,48,32},{31,47,29},{17,36,30},{16,37,35},{34,51,36},{35,50,33},{19,39,34},{21,41,39},{38,53,37},{24,43,41},{40,54,38},{23,44,43},{42,55,40},{28,46,42},{27,47,46},{45,57,44},{32,49,45},{31,50,49},{48,58,47},{36,52,48},{35,53,52},{51,59,50},{39,54,51},{41,56,53},{43,57,56},{55,59,54},{46,58,55},{49,59,57},{52,56,58}};
                BQ.insert(G,I, LaunchCtx(), LaunchPolicy::SYNC, false);
                ID_queue.pop();
            }
            ++I;
            if(I % 10000 == 0){
                progress_bar.update_progress((float)I/(float)num_fullerenes.find(N)->second, "F: " + to_string(I) + "  S: " + to_string(BQ.get_size()));
            }
                
            if(ID_queue.empty()) break;
        }

        BQ.refill_batch(d_validation);
        isomerspace_tutte::tutte_layout(d_validation);
        isomerspace_X0::zero_order_geometry(d_validation, 4.0f);
        cuda_io::reset_convergence_statuses(d_validation);
        for (int i = 0; i < N*10; i+=10){
            isomerspace_forcefield::optimize_batch<BUSTER>(d_validation, 10, N*10);
        }
        
        CuArray<device_real_t> bond_rms(batch_size);
        CuArray<device_real_t> angle_rms(batch_size);
        CuArray<device_real_t> dihedral_rms(batch_size);
        CuArray<device_real_t> bond_max(batch_size);
        CuArray<device_real_t> angle_max(batch_size);
        CuArray<device_real_t> dihedral_max(batch_size);
        CuArray<device_real_t> bond_mean(batch_size);
        CuArray<device_real_t> angle_mean(batch_size);
        CuArray<device_real_t> dihedral_mean(batch_size);
        CuArray<device_real_t> energies(batch_size);
        CuArray<device_real_t> gradient_norm(batch_size);
        CuArray<device_real_t> gradient_rms(batch_size);
        CuArray<device_real_t> gradient_max(batch_size);

        cuda_io::copy(h_validation, d_validation);
        cuda_io::sort(h_validation);
        cuda_io::copy(d_validation, h_validation);
        isomerspace_forcefield::get_bond_rms<BUSTER>(d_validation, bond_rms);
        isomerspace_forcefield::get_angle_rms<BUSTER>(d_validation, angle_rms);
        isomerspace_forcefield::get_dihedral_rms<BUSTER>(d_validation, dihedral_rms);
        isomerspace_forcefield::get_bond_max<BUSTER>(d_validation, bond_max);
        isomerspace_forcefield::get_angle_max<BUSTER>(d_validation, angle_max);
        isomerspace_forcefield::get_dihedral_max<BUSTER>(d_validation, dihedral_max);
        isomerspace_forcefield::get_energies<BUSTER>(d_validation, energies);
        isomerspace_forcefield::get_gradient_norm<BUSTER>(d_validation, gradient_norm);
        isomerspace_forcefield::get_gradient_rms<BUSTER>(d_validation, gradient_rms);
        isomerspace_forcefield::get_gradient_max<BUSTER>(d_validation, gradient_max);

        {   
            
            ofstream file1(path + to_string(N) + "_bond_mean.bin", std::ios::binary);
            ofstream file2(path + to_string(N) + "_bond_rms.bin", std::ios::binary);
            ofstream file3(path + to_string(N) + "_bond_max.bin", std::ios::binary);
            ofstream file4(path + to_string(N) + "_angle_mean.bin", std::ios::binary);
            ofstream file5(path + to_string(N) + "_angle_rms.bin", std::ios::binary);
            ofstream file6(path + to_string(N) + "_angle_max.bin", std::ios::binary);
            ofstream file7(path + to_string(N) + "_dihedral_mean.bin", std::ios::binary);
            ofstream file8(path + to_string(N) + "_dihedral_rms.bin", std::ios::binary);
            ofstream file9(path + to_string(N) + "_dihedral_max.bin", std::ios::binary);
            ofstream file10(path + to_string(N) + "_gradient_norm.bin", std::ios::binary);
            ofstream file11(path + to_string(N) + "_gradient_rms.bin", std::ios::binary);
            ofstream file12(path + to_string(N) + "_gradient_max.bin", std::ios::binary);
            ofstream file13(path + to_string(N) + "_energies.bin", std::ios::binary);
            ofstream file14(path + to_string(N) + "_iterations.bin", std::ios::binary);

            file1.write((char*)bond_mean.data, bond_mean.size_*sizeof(float));
            file2.write((char*)bond_rms.data, bond_rms.size_*sizeof(float));
            file3.write((char*)bond_max.data, bond_max.size_*sizeof(float));
            
            file4.write((char*)angle_mean.data, angle_mean.size_*sizeof(float));
            file5.write((char*)angle_rms.data, angle_rms.size_*sizeof(float));
            file6.write((char*)angle_max.data, angle_max.size_*sizeof(float));

            file7.write((char*)dihedral_mean.data, dihedral_mean.size_*sizeof(float));
            file8.write((char*)dihedral_rms.data, dihedral_rms.size_*sizeof(float));
            file9.write((char*)dihedral_max.data, dihedral_max.size_*sizeof(float));

            file10.write((char*)gradient_norm.data, gradient_norm.size_*sizeof(float));
            file11.write((char*)gradient_rms.data, gradient_rms.size_*sizeof(float));
            file12.write((char*)gradient_max.data, gradient_max.size_*sizeof(float));
            file13.write((char*)energies.data, energies.size_*sizeof(float));
            file14.write((char*)h_validation.iterations, batch_size*sizeof(size_t));

        }
    }
}