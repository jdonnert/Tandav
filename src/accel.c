/* 
 * Collect all accelerations on particle ipart 
 */

#include "globals.h"
#include "accel.h"
#include "timestep.h"
#include "Gravity/gravity.h"
#include "memory.h"
static void Accel_Gravity();

void Compute_Acceleration()
{
	Profile("Accelerations");

#ifdef GRAVITY
	Accel_Gravity();
#endif // GRAVITY

	Profile("Accelerations");

	return ;
}

static void Accel_Gravity()
{

#ifdef GRAVITY_SIMPLE
	Accel_Gravity_Simple();
#endif

#ifdef GRAVITY_TREE
	Build_Gravity_Tree();

Print_Memory_Usage();
	
	Gravity_Tree_Acceleration();
exit(0);
#endif


	return ;
}

