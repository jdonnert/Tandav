struct bunch_information { 
	Float Pos[3];
	int Npart[NPARTYPE];
	int Next_bunch; 
	int Target_task;
	peanoKey Key;
} *Bunchlist; // These will also be the leafs of the dynamic top node tree

void Domain_Decomposition();
void Init_Domain_Decomposition();

struct Domain_Properties {  
	double Size;		// size of smallest cubic box containing all parts
	double Origin[3];	// origin of smallest cubic box containing all parts
} Domain;

