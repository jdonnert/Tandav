#include "io.h"

void Read_ICs(int argc, char *argv[])
{
	Profile("Read");

	switch (Param.Start_Flag) {

	case READ_IC:

		Read_Snapshot(Param.Input_File); // also init particle structures

		break;

	case READ_RESTART:

		Read_Restart_File();

		break;

	case READ_SNAP:

		Assert(argc > 3, "Missing snapshot number in program invokation");    

		Restart.Snap_Counter =  atoi(argv[3]);

		char snap_file[CHARBUFSIZE] = {""};

		sprintf(snap_file, "%s_%03d", Param.Output_File_Base, atoi(argv[3]));

		Read_Snapshot(snap_file);

		break;

	default:

		Assert(false, "Start Flag %d not handled", Param.Start_Flag);

		break;
	}

	Profile("Read");

	return ;
}

unsigned int Largest_Block_Member_Nbytes()
{
	int imax = 0;
	size_t nBytes_max = 0;

	for (int i = 0; i < NBlocks; i++) {
	
		size_t nBytes = Block[i].Ncomp * Block[i].Nbytes;

		if (nBytes > nBytes_max) {

			imax = i;
			nBytes_max = nBytes;
		}
	}

	return Block[imax].Ncomp * Block[imax].Nbytes;
}

/* 
 * Finds the number of particles for block i 
 */

unsigned int Npart_In_Block(const int i, const int *nPart)
{
	unsigned int list[NPARTYPE + 1] = { 0 };

	for (int i = 0; i < NPARTYPE; i++) {

		list[0] += nPart[i];

		list[i+1] = nPart[i];
	}

	int j = (int) Block[i].Target;

	return list[j];
}
