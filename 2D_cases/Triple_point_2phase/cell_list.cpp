#include "glbfunc.h"
#include "particle.h"

#include "level_infor.h"
#include "cell_list.h"


/***************************************************/
/*                                                 */
/*     Functions defined in class "Cell_list"      */
/*                                                 */
/***************************************************/

//-------------------------------------------------------
// initialze all the necessary parameters for cell_list
//-------------------------------------------------------
void Cell_list::Initialize
(Level_info *level_info)
{
  level = level_info->level;
  particle_list.clear();
}
//-------------------------------------------------------
// reset the counter to 0
//-------------------------------------------------------
void Cell_list::Reset_tags()
{
  particle_list.clear();

}
//-------------------------------------------------------
// add particle to cell_list
//-------------------------------------------------------
void Cell_list::Add_particle (Particle *current_particle)
{
  particle_list.push_back(current_particle);
}

//-------------------------------------------------------
// Remove_particle //2026.9.21
//-------------------------------------------------------
void Cell_list::Remove_particle(p_Particle cp)
{
    concurrent_vector<p_Particle> new_list;
    new_list.clear();

    for (int i = 0; i < int(particle_list.size()); ++i)
    {
        if (particle_list[i] != cp)
            new_list.push_back(particle_list[i]);
    }

    particle_list.swap(new_list);
}