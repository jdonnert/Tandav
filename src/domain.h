struct Bunch_Info { 
	shortKey Key;
	float Pos[3];
	int Target; 		// target >=0 -> MPI Rank, target<0 -> -(ipart+1)
	int Npart;
	float CPU_Cost;
} *B; // These will also be the leafs of the dynamic top node tree

//int *Tree2Bunch = NULL, *Bunch2Tree = NULL;

struct Domain_Properties {  
	double Size;		// size of smallest cubic box containing all parts
	double Origin[3];	// origin of smallest cubic box containing all parts
	double Center[3];	// center of smallest cubic box containing all parts
} Domain;

int NBunches;
int NLocal_Bunches;

void Domain_Decomposition();
void Init_Domain_Decomposition();



