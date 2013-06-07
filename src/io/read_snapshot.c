#include "../globals.h"
#include "io.h"

static size_t safe_fread(void *, size_t size, size_t nWanted, FILE *stream);
static void read_file(char *filename, int ReadTask, int LastTask);
static int find_files(char *filename);
static bool find_endianess(FILE * fp);
static int find_block(FILE * fp, const char label[4]);
static void read_header_data(FILE *fp);
static void read_file(char *, const int, const int, const bool swapEndian);

#define SKIP_FORTRAN_RECORD safe_fread(&Fortran_Record, 4, 1, fp, swapEndian);
static int32_t Fortran_Record; // holds the 4 byte Fortran data

void Read_Snapshot(char *input_name)
{
	const nTask = Task.NTask;
	int nIOTasks = Param.No_Output_Files;

	int restFiles = 0;
 	bool swapEndian = false;
	char filename[CHARBUFSIZE];

	rprintf("\nParallel Reading on %d tasks\n",	nIOTasks);

	if (!Task.Rank) {

		restFiles = find_files(input_name);

		strncpy(filename, input_name, CHARBUFSIZE);

		if (restFiles > 1)
        	sprintf(filename, "%s.0", input_name);

		FILE *fp = fopen(filename, "r");

		swapEndian = find_endianess(fp);

		read_header_data(fp, swapEndian); // fills global variable Sim

		fclose(fp);
	}

	MPI_Bcast(&restFiles, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&swapEndian, sizeof(swapEndian), MPI_BYTE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&Sim, sizeof(Sim), MPI_BYTE, 0, MPI_COMM_WORLD);

	while (restFiles > 0) {
		
		if (restFiles >= nTask) { // read in nIO blocks, no communication
		
			int fileNo = Task.Rank + (restFiles - nTask); // backwards

			sprintf(filename, "%s.i", input_name, fileNo);

			for (int i = 0; i < groupSize; i++) {

				if (Task.Rank == groupMaster + i)
					read_file(filename, Task.Rank, Task.Rank, swapEndian);

				MPI_Barrier(MPI_COMM_WORLD);
			}

			restFiles -= nTask;
		} else {
		
			if (nIOTasks > restFiles)
				nIOTasks = restFiles;
			
			int groupSize = nTask / nIOTasks;

			if (nTask % nIOTasks) 
				groupSize++;

			int groupMaster = (Task.Rank / groupSize) * groupSize;

			int groupLast = groupMaster + groupSize - 1;

			if (groupLast > nTask - 1) 
				groupLast = nTask - 1;
		 
			int fileNo = groupMaster / groupSize + restFiles - nIOTasks;

			sprintf(filename, "%s.%i", input_name, fileNo);

			read_file(filename, groupMaster, groupLast, swapEndian);

			restFiles -= nIOTasks; 
		}

		MPI_Barrier(MPI_COMM_WORLD);
	}

	rprintf("Reading completed\n");
	

// find snapshot type & endianess
// read first file header and fill global vars
// do parallel decomposition
// read :
// 		create communicator
//		<0> header
//		communicate Npart
//		decompose particles
//		make space
//		read data on <0>
//		do one to all communication
//		sort in locally
//		deallocate
	return ;
}


static void read_file(char *filename, const int readTask, const int lastTask, 
		const bool swapEndian)
{
	FILE *fp = NULL;
	int32_t nPart[NO_PART_TYPES] = { 0 }, nTot = 0;
	
	MPI_Comm mpi_comm_read = Create_MPI_Communicator(readTask, lastTask);
	
	if (Task.Rank == readTask) { // get nPart

		fp = fopen(filename, "r");
	 
		find_block(fp, "HEAD");

		SKIP_FORTRAN_RECORD;

		safe_fread(nPart, NO_PART_TYPES*sizeof(*nPart), 6, fp);
		
		SKIP_FORTRAN_RECORD;

		for (int i=0; i<NO_PART_TYPES; i++)
			ntot += nPart[i];
	}

	MPI_Bcast(nPart, NO_PART_TYPES, MPI_FLOAT, 0, mpi_comm_read);
	MPI_Bcast(nTot, 1, MPI_FLOAT, 0, mpi_comm_read);

	rprintf("\nReading file <%s> on Task <%i-%i> \n"
            "   Sph   <%9li>  DM     <%9li>    \n"
            "   Disk  <%9li>  Bulge  <%9li>    \n"
            "   Star  <%9li>  Bndry  <%9li>    \n"
            "   Total <%9li> \n\n",
			filename, ReadTask, LastTask,
			nPart[0], nPart[1], nPart[2],
			nPart[3], nPart[4], nPart[5], 
            ntot); fflush(stdout);
	 

	/* find npart per CPU and type */
	int groupSize = LastTask - ReadTask;

	int32_t nPartGet[NO_PART_TYPES] = { 0 };
	
	int32_t nPartGetTotal = 0;
	
	for (int i = 0; i < N_part_types; i++) {
	
		for (int j=ThisTask.Rank-ReadTask; j<nPart[i]; j+=groupSize) {
		
			nPartGet[i]++;
			
			nPartGetTotal++;
		}
	}


	fclose(fp);
	
	MPI_Comm_free(mpi_comm_read);

	return ;
}
static void read_header_data(FILE *fp, bool swapEndian) 
{ 
	const int ntype = NO_PART_TYPES;
	
	struct gadget_header head;

	rewind(fp);

	int32_t blocksize = find_block(fp, "HEAD");

	Assert(blocksize == 256, "Header corrupted");

	SKIP_FORTRAN_RECORD;

	safe_fread(head.Npart, ntype*sizeof(*head.Npart), 6, fp);
	safe_fread(&head.Massarr, ntype*sizeof(head.Massarr), 6, fp);
	safe_fread(&head.Time, ntype*sizeof(head.Time), 1, fp);
	safe_fread(&head.Redshift, ntype*sizeof(head.Redshift), 1, fp);
	safe_fread(&head.FlagSfr, ntype*sizeof(head.FlagSfr), 1, fp);
	safe_fread(&head.FlagFeedback, ntype*sizeof(head.FlagFeedback), 1, fp);
	safe_fread(head.Nall, ntype*sizeof(*head.Nall), 1, fp);
	safe_fread(&head.FlagCooling, ntype*sizeof(head.FlagCooling), 1, fp);
	safe_fread(&head.NumFiles, ntype*sizeof(head.NumFiles), 6, fp);
	safe_fread(&head.Boxsize, ntype*sizeof(head.Boxsize), 1, fp);
	safe_fread(&head.Omega0, ntype*sizeof(head.Omega0), 1, fp);
	safe_fread(&head.OmegaLambda, ntype*sizeof(head.OmegaLambda), 1, fp);
	safe_fread(&head.HubbleParam, ntype*sizeof(head.HubbleParam), 1, fp);
	safe_fread(&head.FlagAge, ntype*sizeof(head.FlagAge), 1, fp);
	safe_fread(&head.FlagMetals, ntype*sizeof(head.FlagMetals), 1, fp);
	safe_fread(&head.NallHighWord, ntype*sizeof(head.NallHighWord), 6, fp);

	for (i=Sim.NpartTotal=0; i<NO_PART_TYPES; i++) {
		Sim.Mpart[i] = head.Massarr[i];

		Sim.Npart[i] = (uint64_t)head.Npart[i];
		Sim.Npart[i] += ((uint64_t)head.NallHW[i]) << 32;

		Sim.NpartTotal += Sim.Npart[i];
	}

	Sim.Boxsize = head.boxsize; // rest defined in parameter file

	rewind(fp);

	return ;
}

static int find_block(FILE * fp, const char label[4])
{
	int32_t blocksize;

	rewind(fp);

	for (;;) {
		char blocklabel[4] = {"    "};
		
		SKIP_FORTRAN_RECORD;
		
		int nRead = safe_fread(blocklabel, 4 * sizeof(char), 1, fp);

		Assert(nRead == 1, "Format 2 File broken");

		safe_fread(blocksize, 4 * sizeof(char), 1, fp);

		Assert(nRead == 1, "Format 2 File broken");

        SKIP_FORTRAN_RECORD;
			
        if (strncmp(label, blocklabel, 4) == 0) 
			break; // found it
			
		fseek(fp, blocksize, 1); // skip to end of block
		
		blocksize = 8;

		if (feof(fp)) 
			break; // not in this file
	}

	return blocksize-8; // remove 8 bytes from Fortran header
}


static bool find_endianess(FILE *fp) 
{
	rewind(fp);

	int32_t testRecord = 0;
	
	size_t nRead = fread(&testRecord, 4, 1, fp);

	Assert(nRead == 1, "Couldn't read Fortran record");

	bool swapEndian = false;

	if (testRecord == 134217728) { // first block is always 8 bytes in format 2

		printf("\nEnabling Endian Swapping\n");

		swapEndian = true;
	}

	rewind(fp);

	SKIP_FORTRAN_RECORD

	Assert(Fortran_Record == 8, "Binary Fortran File Format Broken");

	fclose(fp);

	return swapEndian;
}

static int find_files(char *filename)
{
	char buf[CHARBUFSIZE] = "\0";

	FILE *fp = fopen(filename, "r");

	int nFiles = 0;

	if (fp != NULL) { 

		nFiles++;

		fclose(fp);

	} else {

		for (;;) {

            sprintf(buf, "%s.%i", filename, nFiles);
			
            if (!(fp = fopen(buf, "r"))) 
				break;
			
			fclose(fp);
			
            nFiles++;

            Assert(nFiles<10000, "Found 10000 files, holy cow !");
		}

    	Assert(nFiles, "Can't open input file as '%s' or '%s'", filename, buf);
	}

	return nFiles;
}

/* this fread wrapper checks for eof, corruption and swaps endianess 
 * if swapEndian == true */
static size_t safe_fread(void *data, size_t size, size_t nWanted, FILE *stream, 
		bool swapEndian)
{
	size_t nRead = fread(data, size, nWanted, stream);

	if (feof(stream)) // we catch EOF further up, because of block finding
		nRead = nWanted = 0;
		
	Assert(nRead == nWanted, "Read %zu bytes, but %zu wanted", nRead, nWanted);

	if (swapEndian && nRead) { // swap endianess

		char buf[16];

		for (int j=0; j < nWanted; j++) {
			memcpy(buf, &(( (char *)data )[j * size]), size);
			
			for (int i=0; i < size; i++)
				( (char*)data )[j*nWanted + i] = buf[size - i - 1];
		}
	}

	return nRead;
}
#undef SKIP_FORTRAN_RECORD
