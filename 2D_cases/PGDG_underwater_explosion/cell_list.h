#ifndef CELL_LIST_H
#define CELL_LIST_H
#include <cmath>
#include "glbcls.h"
#include "Mypool.h"

using namespace std;


class Cell_list{
  int      level;

public:

// variables
  concurrent_vector <p_Particle> particle_list;
  Cell_list(){};

// functions
  void Initialize (Level_info *level_info);
  void Reset_tags ();
  void Add_particle (Particle *current_particle);
  void Remove_particle(p_Particle cp);//2026.9.21
};

#endif
