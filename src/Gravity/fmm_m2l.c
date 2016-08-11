#include "fmm.h"

#ifdef GRAVITY_FMM

static void prepare_M2L();
static void push(const int src, const int dst);
static void pop(int*, int*);
static void interact_FMM(const int, const int);
static bool MAC(const int a, const int b);
static void p2p_kernel(const int a, const int b);
static void m2l_kernel(const int a, const int b);

static struct Work_Stack {
	int src;
	int dst;
	int src_rank;
	int dst_rank;
} *wstack;

static int N = 0;

static omp_lock_t wstack_lock;

/* Dual tree walk FMM (Dehnen 2001, Yokota 2012) */

void Gravity_FMM_M2L()
{
	Profile("FMM_M2L");

	prepare_M2L();

	push(0,0);

	while (N > 0) { // Yokata 2012 Alg 3
		
		int src = 0, dst = 0;

		pop(&src, &dst);

		int parent = dst;

		if (Level(FMM, src) > Level(FMM, dst))
			parent = src;

		int run = parent + 1; // leaf is never pushed
		int sibling = parent + FMM.DNext[dst];
			
		while (run < sibling) { // loop over children
	
			#pragma omp task
			interact_FMM(src, run);

			run += imax(1, FMM.DNext[run]);
		}
	
	} // while N > 0

	Free(wstack);

	Profile("FMM_M2L");

	return ;
}

void M2L_Setup()
{
	omp_init_lock(&wstack_lock);

	return ;
}

static void prepare_M2L()
{
	size_t nBytes = Task.Npart_Total_Max*sizeof(*wstack);
	
	wstack = Malloc(nBytes, "wstack");

	memset(wstack, 0, nBytes);

	N = 0;
	
	return ;
}

static void push(const int src, const int dst)
{
	omp_set_lock(&wstack_lock);

	wstack[N].src = src;
	wstack[N].dst = dst;

	N++;

	omp_unset_lock(&wstack_lock);

	return ;
}

static void pop(int *src, int *dst)
{
	omp_set_lock(&wstack_lock);

	*src = wstack[N].src;
	*dst = wstack[N].dst;

	wstack[N].src = wstack[N].dst = 0;

	N--;

	omp_unset_lock(&wstack_lock);

	return ;
}

/* Yokata 2012, Alg 4 */

static void interact_FMM(const int a, const int b)
{
	if ((FMM.DNext[a] < 0) && (FMM.DNext[b] < 0))
		p2p_kernel(a,b);
	else if (MAC(a,b)) 
		m2l_kernel(a,b);
	else
		push(a,b);

	return;
}

/* Multipole Acceptance Criterion (Dehnen 2001) */

static bool MAC(const int a, const int b)
{
	Float dr[3] = { FMM.CoM[0][a] - FMM.CoM[0][b],
					FMM.CoM[1][a] - FMM.CoM[1][b],
					FMM.CoM[2][a] - FMM.CoM[2][b]};

	//Periodic_Nearest(dr); // PERIODIC 
	
	Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];

	bool result = false;

	if (r2 > p2(FMM.Rcrit[a] + FMM.Rcrit[b])) // (Dehnen 2002, eq 10) 
		result = true;

	return result;
}

static void p2p_kernel(const int a, const int b)
{

	return ;
}

static void m2l_kernel(const int a, const int b)
{
	Float dr[3] = { FMM.CoM[0][a] - FMM.CoM[0][b],
					FMM.CoM[1][a] - FMM.CoM[1][b],
					FMM.CoM[2][a] - FMM.CoM[2][b] } ;
		
	//Periodic_Nearest(dr); // PERIODIC 
	
	Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]; // assume r > epsilon
	
	Float force[3] = { (FMM.Mass[a] + FMM.Mass[b]) * dr[0]/r2,
					   (FMM.Mass[a] + FMM.Mass[b]) * dr[1]/r2,
					   (FMM.Mass[a] + FMM.Mass[b]) * dr[2]/r2 };

	FMM.Force[0][a] = force[0];
	FMM.Force[1][a] = force[1];
	FMM.Force[2][a] = force[2];

	FMM.Force[0][b] = -force[0];
	FMM.Force[1][b] = -force[1];
	FMM.Force[2][b] = -force[2];

	return ;
}
#endif // GRAVITY_FMM

// Copyright (C) 2013 Julius Donnert (donnert@ira.inaf.it)
