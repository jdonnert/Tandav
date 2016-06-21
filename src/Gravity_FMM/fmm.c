#include "fmm.h"

#ifdef GRAVITY_FMM

void Gravity_Acceleration() 
{
	//if (Sig.Tree_Update)
		Gravity_FMM_Build(); 

	Gravity_FMM_P2L();

	//Gravity_FMM_M2L();
	
	//Gravity_FMM_L2L();
	
	//Gravity_FMM_P2P();

	return ;
}


#endif // GRAVITY_FMM
