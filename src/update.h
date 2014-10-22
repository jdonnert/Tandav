enum Update_Parameters {
	BEFORE_MAIN_LOOP,
	BEFORE_STEP,
	BEFORE_FIRST_KICK,
	BEFORE_SNAPSHOT,
	BEFORE_DRIFT,
	FORCES,
	BEFORE_SECOND_KICK,
};

void Update(enum Update_Parameters stage); 
