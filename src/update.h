void Update(enum Update_Parameters stage); 

enum Update_Parameters {
	BEFORE_MAIN_LOOP,
	AFTER_FIRST_KICK,
	AFTER_DRIFT,
	AFTER_SECOND_KICK
};
