static struct bunch_information { 
	Float Pos[3];
	int Npart[NPARTYPE];
	int Next_bunch; 
	int Target_task;
	peanoKey Key;
} *Bunchlist; // These will also be the leafs of the dynamic top node tree

void Domain_Decomposition();
void Init_Domain_Decomposition();

double Domain_Corner[3];
