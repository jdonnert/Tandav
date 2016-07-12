/*
 * Automatic particle handling. Because of the changed data structures, we
 * can't simply allocate, shrink and enlarge the particles. Instead
 * we now have to loop over all its members and their components and change
 * all these arrays individually. The necessary information is in a const
 * struct array in particles_fields.h, which we may autogenerate from 
 * includes.h at some point
 */

#include "particles.h"

struct Particle_Data P = { NULL };
struct Gas_Particle_Data G = { NULL };

size_t sizeof_P = 0; 
const int NP_Fields = ARRAY_SIZE(P_Fields);

static void find_particle_sizes();

static omp_lock_t Particle_Lock; // the big bad particle lock

/*
 * We loop through the particle structures with two running pointers. One 
 * points through the structure (pointers to adresses are all the same 
 * size: 64 bit), the other points through the memory allocated.
 */

void Allocate_Particle_Structures()
{
	#pragma omp parallel // Task is threadprivate
	{

	const double npart_per_rank = (double) Sim.Npart_Total/(double) NRank;

	Task.Npart_Total_Max = ceil(npart_per_rank * Param.Part_Alloc_Factor);

	for (int i = 0; i < NPARTYPE; i++)
		Task.Npart_Max[i] = (double) Sim.Npart[i]/NRank 
							* Param.Part_Alloc_Factor;

	} // omp parallel

	find_particle_sizes(); 

	size_t nBytes = Task.Npart_Total_Max * sizeof_P;

	rprintf("\nReserving space for %llu particles per task in *P,"
			" factor %g\n", Task.Npart_Total_Max, Param.Part_Alloc_Factor);

	omp_init_lock(&Particle_Lock);

	omp_set_lock(&Particle_Lock);

	void * restrict * run_P = (void * restrict *) &P.Type; // first field in P

	for (int i = 0; i < NP_Fields; i++) {
			
		int nComp = P_Fields[i].N;

		nBytes = Task.Npart_Total_Max * P_Fields[i].Bytes;
		
		for (int j = 0; j < nComp;  j++) {

			char name[CHARBUFSIZE] = { "" }; 

			sprintf(name, "P.%s[%d]", P_Fields[i].Name, j);

			*run_P = Malloc(nBytes, name); // memory block

			run_P++; // next field in P is 8 bytes or one pointer away
		}
	}
	
	//G = Malloc(Task.Npart_Max[0] * sizeof(*G), "G");

	omp_unset_lock(&Particle_Lock);

	Print_Memory_Usage();

	return ;
}

/*
 * This used to be easier
 */

static void find_particle_sizes()
{
	sizeof_P = 0;

	for (int i = 0; i < NP_Fields; i++) {

		sizeof_P += P_Fields[i].Bytes * P_Fields[i].N;
		//sizeof_G += G_Fields[i].Bytes * G_Fields[i].N;
		//....
	}
	
	printf("\nsizeof(P) = %zu byte\n", sizeof_P);

	return ;
}

/* 
 * Reallocate the Particle structure. Takes the relative change
 * as argument, not the total number. Add or Remove via sign of nPart.
 * Also updates Task.Npart and Task.NPartTotal. 
 * Expands P so that space for nPart[type] is at offset[type]
 * Contracts P so that the last nPart[type] particles are removed 
 * Note that actually no real allocation is taking place, because
 * that would fragment memory. Instead this needs to stay with in the
 * limit set by Param.Part_Alloc_Factor. 
 */

void Reallocate_P_Info(const char *func, const char *file, int line,
		const int dNpart[NPARTYPE], size_t offset_out[NPARTYPE])
{
	
	#pragma omp single
	for (int i = 0; i < NPARTYPE; i++)
		Assert(Task.Npart[i] + dNpart[i] <= Task.Npart_Max[i],
			"Too many particles type %d on this task. \n"
			"Have %d, want %d, max %d \nCurrent Part alloc factor = %g",
			i, Task.Npart[i], dNpart[i], Task.Npart_Max[i], 
			Param.Part_Alloc_Factor);

	int offset[NPARTYPE] = { 0 }, new_npart_total = 0;
	int new_npart[NPARTYPE] = { 0 };

	#pragma omp single
	for (int type = 0; type < NPARTYPE; type++) { // calc offset

		new_npart[type] = Task.Npart[type] + dNpart[type];

		new_npart_total += new_npart[type];

        Assert(new_npart[type] >= 0, "Can't alloc negative particles,"
			" type %d, delta %d, current %d,\n"
			"requested from %s, %s(), line %d",
			type, dNpart[type], Task.Npart[type], file, func, line);

		if (dNpart[type] == 0)
			continue; // don't need offset here

		for (int i = 0; i <= type; i++)
			offset[type] += new_npart[i];

		offset[type] -= MAX(0, dNpart[type]); // correct for dNpart>0
	}

	int nMove = Task.Npart_Total; // move particles left

	#pragma omp single
	for (int type = 0; type < NPARTYPE; type++) {

		nMove -= Task.Npart[type];

		if (dNpart[type] >= 0 || Task.Npart[type] == 0 || nMove == 0)
			continue;

		int src = offset[type] + fabs(dNpart[type]);
		int dest = offset[type];

		for (int i = 0; i < NP_Fields; i++) {
		
			for (int j = 0; j < P_Fields[i].N; j++) {

				char *beg = Select_Particle(i, j, src);
				char *end = Select_Particle(i, j, dest);

				memmove(end, beg, nMove*P_Fields[i].Bytes);
			} // j
		} // i
	}

	nMove = Task.Npart_Total; // move particles right

	#pragma omp single
	for (int type = 0; type < NPARTYPE-1; type++) {

		nMove -= Task.Npart[type];

		if (dNpart[type] <= 0 || Task.Npart[type] == 0 || nMove == 0)
			continue;

		int src = offset[type];
		int dest = offset[type] + dNpart[type];
	
		for (int i = 0; i < NP_Fields; i++) {
		
			for (int j = 0; j < P_Fields[i].N; j++) {

				char *beg = Select_Particle(i, j, src);
				char *end = Select_Particle(i, j, dest);

				memmove(end, beg, nMove*P_Fields[i].Bytes);
			}
		}
	}

	#pragma omp parallel  // book-keeping
	{

	Task.Npart_Total = new_npart_total;

	for (int type = 0; type < NPARTYPE; type++)
		Task.Npart[type] = new_npart[type];

	} // parallel 

	if (offset_out != NULL) // return ptrs to freed space
		memcpy(offset_out, offset, NPARTYPE*sizeof(*offset));

	return ;
}

/*
 * Returns a pointer to particle "ipart", field "field", component "comp".
 * This lets us move particles in an automated way.
 * Pointer fun for the whole family.
 */

char * Select_Particle(const size_t field, const int comp, const int ipart)
{
	Assert(comp >= 0 && comp < P_Fields[field].N, 
			"Component %d not valid in field %d", comp, field);

	Assert(field < NP_Fields, "Field %d does not exist", field);

	Assert(ipart >= 0 && ipart < Task.Npart_Total, 
			"ipart %d out of range", ipart);

	char * restrict * result = (char * restrict *) &P.Type;

	for (int i = 0; i < field; i++)
		for (int j = 0; j < P_Fields[i].N; j++)
			result++;
	
	result += comp; // now points to P.$field[$comp]

	return *result + ipart*P_Fields[field].Bytes;
}

Copyright (C) 2013 Julius Donnert (donnert@ira.inaf.it)
