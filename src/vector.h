void Find_Vectors();
void Setup_Vectors();

extern struct Vector_Data {
	int * restrict First;
	int * restrict Last;
	peanoKey * restrict Key;
} Vec;

int NVec;
