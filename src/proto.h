/* Dump all global (!) function prototypes here 
 * I.e. functions down to main + 1 */
enum Update_Parameters {
	BEFORE_MAIN_LOOP,
	AFTER_FIRST_KICK,
	AFTER_DRIFT,
	AFTER_SECOND_KICK
};

extern void Read_Parameter_File(char *);
extern void Write_Parameter_File(char *);
extern void Setup();
extern void Read_Snapshot();
extern void Read_Restart_File();
extern void Init();
extern void Update(enum Update_Parameters);
extern void Kick_First_Halfstep();
extern void Drift();
extern void Write_Snapshot();
extern void Write_Restart_File();
extern void Kick_Second_Halfstep();

extern void print_compile_time_settings();
extern void Qsort(void const *, size_t, size_t, 
		int (*)(const void*, const void*));


