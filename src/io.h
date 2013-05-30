struct parameter_definitions {
	char tag[MAXSTRINGLENGTH];
	char val[MAXSTRINGLENGTH];
	void *addr;
	enum param_type {
		FLOAT,
		INT,
		STRING,
	} id;
};

static struct parameter_definitions pDef[] = { 
	/* tag, 	std val, address, 		 type */
	{"Boxsize", "10000", &Param.Boxsize, FLOAT},
	{"No_IOTasks", "1", &Param.No_Output_Files, INT},
	{"Input_File", "IC_file", &Param.Input_File, STRING},
	{"Output_File_Base", "snap_", &Param.Output_File_Base, STRING}
};

static const int nTags = sizeof(pDef) / sizeof(*pDef);
