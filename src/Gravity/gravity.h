
#if defined(GRAVITY) && defined(GRAVITY_SIMPLE)
void Gravity_Simple_Accel();
#else
inline void Gravity_Simple_Accel() {};
#endif // GRAVITY && GRAVITY_SIMPLE

#if defined(GRAVITY) && defined(GRAVITY_TREE)

void Setup_Gravity_Tree();
void Gravity_Tree_Build();
void Gravity_Tree_Acceleration();
void Gravity_Tree_Update_Kicks(const int ipart, const double dt);
void Gravity_Tree_Update_Topnode_Kicks();
void Gravity_Tree_Update_Drift(const double dt);

extern struct Tree_Node {
	int DNext;			// Distance to the next node; or particle -DNext-1
	uint32_t Bitfield;	// bit 0-5:level, 6-8:key, 9:local, 10:top, 11-31:free
	int DUp;			// Number of nodes to the parent
	int Npart;			// Number of particles in node
	Float Pos[3];		// Node Center
	Float Mass;			// Total Mass of particles inside node
	Float CoM[3];		// Center of Mass
	Float Dp[3];		// Velocity of Center of Mass
} * restrict Tree;

double Epsilon[NPARTYPE], Epsilon2[NPARTYPE]; // grav. softening

#else

inline void Setup_Gravity_Tree() {}; 
inline void Gravity_Tree_Build() {};
inline void Gravity_Tree_Acceleration() {};
inline void Gravity_Tree_Update_Kicks(const int ipart, const double dt) {};
inline void Gravity_Tree_Update_Topnode_Kicks() {};
inline void Gravity_Tree_Update_Drift(const double dt) {};

#endif // ! (GRAVITY && GRAVITY_TREE)

#if defined(GRAVITY) && defined(GRAVITY_FMM)

inline void Setup_Gravity_FMM();
inline void Free_Gravity_FMM();
inline void Gravity_FMM_Build();
inline void Gravity_FMM_Acceleration();
inline void Gravity_FMM_Update_Kicks();
inline void Gravity_FMM_Update_Topnode_Kicks();
inline void Gravity_FMM_Update_Drift(const double dt);

#else

inline void Setup_Gravity_FMM() {};
inline void Free_Gravity_FMM() {};
inline void Gravity_FMM_Build() {};
inline void Gravity_FMM_Acceleration() {};
inline void Gravity_FMM_Update_Kicks() {};
inline void Gravity_FMM_Update_Topnode_Kicks() {};
inline void Gravity_FMM_Update_Drift(const double dt) {};

#endif // GRAVITY && GRAVITY_FMM

