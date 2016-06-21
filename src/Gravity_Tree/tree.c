#include "tree.h"

#ifdef GRAVITY_TREE

/*
 * Barnes & Hutt gravity tree driver routine
 */

void Gravity_Acceleration()
{
	if (Sig.Tree_Update)
		Gravity_Tree_Build();

	if (Sig.Prepare_Step) {

		Gravity_Tree_Walk(true); // with BH criterion
	
		Safe_Last_Accel();
	}

	Gravity_Tree_Walk(false);

	return ;
}

#endif
