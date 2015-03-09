#if defined(GRAVITY) && defined(GRAVITY_TREE)
int Level(const int node); // bitfield functions

enum Tree_Bitfield { LOCAL=9, TOP=10, UPDATED=11 }; // offset by one

Float Node_Size(const int node);
bool Node_Is(const enum Tree_Bitfield bit, const int node);
void Node_Set(const enum Tree_Bitfield bit, const int node);
void Node_Clear(const enum Tree_Bitfield bit, const int node);
#endif
