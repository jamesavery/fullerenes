//#include <cstdlib>
#include <vector>
//#include <math.h>
#include <iostream>

void connect(const int &i, const int &j, int connectivity[][6], std::vector<int> &uv)
{
  connectivity[i][uv[i]] = j+1;
  connectivity[j][uv[j]] = i+1;
  ++uv[i];
  ++uv[j];
}

bool windup_general(const std::vector<int> pot_spiral, const unsigned int size, std::vector<int> &pos, std::vector<int> &dist, int connectivity[][6]){
  //number of used valencies per vertex (0-5(6))
  std::vector<int> used_valencies(size);
  //list of vertices that have open valencies
  std::vector<int> open_valencies;
  // number of used pentagons (0-12)
  int used_pentagons;

  //connect first two faces
  connectivity[0][0] = 2;
  connectivity[1][0] = 1;
  used_valencies[0] = 1;
  used_valencies[1] = 1;
  open_valencies.push_back(0);
  open_valencies.push_back(1);
  used_pentagons = 12 - pot_spiral[0] - pot_spiral[1];
  int x; //jump distance (x=1 is no jump)
  
  //iterate over atoms
  //k=0, k=1 have been done already
  for (int k=2; k<size; ++k){
    used_pentagons += 6 - pot_spiral[k];
  
    //check if jump
    x = 0;
    if(pot_spiral[open_valencies.front()] - used_valencies[open_valencies.front()] == 2 && open_valencies.size() > 6){
      while(pot_spiral[*(open_valencies.begin()+x)] - used_valencies[*(open_valencies.begin()+x)] == 2){
        ++x;
      }
      //two error cases
      if (pot_spiral[*(open_valencies.begin()+x)] - used_valencies[*(open_valencies.begin()+x)] == 1 &&
        pot_spiral[*(open_valencies.begin()+x+1)] - used_valencies[*(open_valencies.begin()+x+1)] == 1 &&
        pot_spiral[*(open_valencies.begin()+x+2)] - used_valencies[*(open_valencies.begin()+x+2)] == 1 &&
        pot_spiral[*(open_valencies.begin()+x+3)] - used_valencies[*(open_valencies.begin()+x+3)] == 1)
                {std::cout << "There is no spiral." << std::endl;
        return 1;}
      if (pot_spiral[*(open_valencies.begin()+x)] - used_valencies[*(open_valencies.begin()+x)] == 1 &&
        pot_spiral[*(open_valencies.begin()+x+1)] - used_valencies[*(open_valencies.begin()+x+1)] == 1 &&
        pot_spiral[*(open_valencies.begin()+x+2)] - used_valencies[*(open_valencies.begin()+x+2)] == 1 &&
        pot_spiral[*(open_valencies.begin()+x)] == 5)
                {std::cout << "There is no spiral." << std::endl;
        return 1;}
      //two jump cases
      if ((pot_spiral[*(open_valencies.begin()+x)] - used_valencies[*(open_valencies.begin()+x)] == 1 &&
         pot_spiral[*(open_valencies.begin()+x+1)] - used_valencies[*(open_valencies.begin()+x+1)] == 1 &&
         pot_spiral[*(open_valencies.begin()+x)] == 5) || 
        (pot_spiral[*(open_valencies.begin()+x)] - used_valencies[*(open_valencies.begin()+x)] == 1 &&
        pot_spiral[*(open_valencies.begin()+x+1)] - used_valencies[*(open_valencies.begin()+x+1)] == 1 &&
        pot_spiral[*(open_valencies.begin()+x+2)] - used_valencies[*(open_valencies.begin()+x+2)] == 1 &&
        pot_spiral[*(open_valencies.begin()+x)] == 6))
        ++x;
        else
        x=0;
    }

    //jump positions and distances
    if(x>1){// x=1 is no jump
      pos.push_back(k);
      dist.push_back(x);
    }

    // perform cyclic rotation on open_valencies
    for(int i = 1; i<x; ++i){
      int j = open_valencies.front();
      open_valencies.erase(open_valencies.begin());
      open_valencies.push_back(j);
    }

    //connect k to k-1
    connect(open_valencies.back(), k, connectivity, used_valencies);

    //connect k to k-2, etc
    while(open_valencies.size() != 0 && pot_spiral[open_valencies.back()] - used_valencies[open_valencies.back()] == 0){
//      std::cout << pot_spiral[open_valencies.back()] << used_valencies[open_valencies.back()] << std:: endl;
      open_valencies.erase(open_valencies.end()-1);
      if(open_valencies.size() != 0 && pot_spiral[k] - used_valencies[k] != 0){
        connect(open_valencies.back(), k, connectivity, used_valencies);
      }
    }

//    std::cout << k << "vii" << std::endl;
    //connect k to oldest unconnected (etc)
    if(pot_spiral[k] - used_valencies[k] != 0){
      connect(open_valencies.front(), k, connectivity, used_valencies);
    }
//    std::cout << open_valencies.front() << " " << pot_spiral[open_valencies.front()] << " " << used_valencies[open_valencies.front()] << std:: endl;
    while(open_valencies.size() != 0 && pot_spiral[open_valencies.front()] - used_valencies[open_valencies.front()] == 0){
      open_valencies.erase(open_valencies.begin());
      if(open_valencies.size() != 0 && pot_spiral[k] - used_valencies[k] != 0){
        connect(open_valencies.front(), k, connectivity, used_valencies);
      }
    }
    
    //append k to the open valencies (after making sure it has some valencies left)
    if (pot_spiral[k] - used_valencies[k] != 0){//the current atom is not saturated (which may only happen for the last one)
      open_valencies.push_back(k);
    }else{
      if(k + 1 != size){
        std::cout << "Fail 1 (cage closed but faces left)" << std::endl;
        return 1;
      }
    }
  }//iterate over atoms
  if(open_valencies.size() != 0){
      std::cout << "Fail 2 (cage not closed but no faces left)" << std::endl;
      return 1;
  } 
  return 0;//success
}//windup_general


int main(int argc, char* argv[])
{
  unsigned int size = 192;
//  unsigned int size = 16;
  std::vector<int> potential_spiral (size,6);
//normal
//  int pentagon_indices[] = {1,2,3,4,5,7,10,12,13,14,15,16};
//jump
//  int pentagon_indices[] = {1,2,3,5,8,9,11,12,13,14,15,16};
//jump, 380
  int pentagon_indices[] = {1,2,3,137,138,147,148,157,158,166,173,180};
  std::cout << "The pentagon indices: ";
  for (int i=0; i<12; ++i)
  {
    int index = pentagon_indices[i];
    potential_spiral[index-1] = 5;
    std::cout << index << ", ";
  }
  std::cout << std::endl << std::endl;

    // print connectivities
  //std::vector<std::vector<int> > connectivities;
  int connectivity[size][6];
    for (int i=0; i < size; ++i){
      for (int j=0; j < 6; ++j){
        connectivity[i][j] = 0;
      }
    }

std::vector<int> jump_positions;
std::vector<int> jump_distances;

if (windup_general(potential_spiral, size, jump_positions, jump_distances, connectivity))
{return 1;}

for (int i=0; i < size; ++i){
  for (int j=0; j < 6; ++j){
    std::cout << connectivity[i][j] << ", " ;
  }
  std::cout << std::endl;
}

std::cout << jump_positions.size() << " jumps required.";
for (std::vector<int>::iterator i(jump_positions.begin()); i<jump_positions.end(); ++i)
{
  std::cout << *i+1 << ", ";//because k is relative to 0
}

for (std::vector<int>::iterator i(jump_distances.begin()); i<jump_distances.end(); ++i)
{
  std::cout << *i << ", ";
}
std::cout << std::endl;

  return 0;
}

