struct bunch_information { 
	int Target_Rank;
	int First_Particle;
	int Npart;
	peanoKey Key;
	int Global_Tree_Offset;
	float CPU_Cost;
	Float CoM[3];
	Float Mass;
} *Bunch; // These will also be the leafs of the dynamic top node tree

int NBunches;
int NLocal_Bunches;

void Domain_Decomposition();
void Init_Domain_Decomposition();

struct Domain_Properties {  
	double Size;		// size of smallest cubic box containing all parts
	double Origin[3];	// origin of smallest cubic box containing all parts
	double Center[3];	// center of smallest cubic box containing all parts
} Domain;

