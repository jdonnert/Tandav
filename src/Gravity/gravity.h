#ifdef GRAVITY

void Accel_Gravity_Simple(const int ipart, double *acc, double *pot);

#ifdef GRAVITY_TREE

void Init_Tree();
void Build_Tree();

#endif // GRAVITY_TREE

#endif // GRAVITY
