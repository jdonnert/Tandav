struct Bunch_Info { // These will also be the leafs of the dynamic top node tree
	shortKey Key;		// Largest Peano key held by this bunch
	float Sample_Pos[3];// Random (!!!) position inside the bunch 
	int Target; 		// target >=0 -> MPI Rank, target<0 -> -(ipart+1)
	int Npart;
	float CPU_Cost;
	int Level;
} *B; 

//int *Tree2Bunch = NULL, *Bunch2Tree = NULL;

struct Domain_Properties {  
	double Size;		// size of smallest cubic box containing all particles
	double Origin[3];	// origin of smallest cubic box containing all particles
	double Center[3];	// center of smallest cubic box containing all partciles
} Domain;

int NBunches;
int NLocal_Bunches;

void Domain_Decomposition();
void Init_Domain_Decomposition();



