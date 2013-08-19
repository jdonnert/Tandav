#include "../globals.h"
#include "io.h"

#define SKIP_FORTRAN_RECORD safe_fread(&Fortran_Record, 4, 1, fp, swapEndian);
static int32_t Fortran_Record; // holds the 4 byte Fortran data

static int safe_fread(void *, size_t, size_t, FILE *, bool);
static int find_files(char *);
static bool find_endianess(FILE *);
static int find_block(FILE *, const char label[4], const bool);
static void read_header_data(FILE *fp, const bool, int);
static void read_file(char *, const bool, const int, const int, MPI_Comm);
static void empty_comm_buffer(char *, const int, const int *, const size_t *);

static void generate_masses_from_header();

static char * restrict RecvBuf = NULL, * restrict ReadBuf = NULL;

void Read_Snapshot(char *input_name)
{
	const int nTask = Sim.NTask;
	
	int nIOTasks = Param.NumIOTasks;

	MPI_Comm mpi_comm_read;

	int nFiles = 0;
	bool swapEndian = false;
	char filename[CHARBUFSIZE] = "";

	if (!Task.Rank) {

		nFiles = find_files(input_name);
		
		rprintf("\nParallel Reading of %d files on %d tasks\n\n", 
				nFiles, nIOTasks);

		strncpy(filename, input_name, CHARBUFSIZE);
		if (nFiles > 1)
        	sprintf(filename, "%s.0", input_name);

		FILE *fp = fopen(filename, "r");

		swapEndian = find_endianess(fp);

		read_header_data(fp, swapEndian, nFiles); // fills global var 'Sim'
		
		fclose(fp);
	}
	
	MPI_Bcast(&nFiles, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&swapEndian, sizeof(swapEndian), MPI_BYTE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&Sim, sizeof(Sim), MPI_BYTE, 0, MPI_COMM_WORLD);
	
	int restFiles = nFiles;

	while (restFiles) {

		nIOTasks = fmin(nIOTasks,restFiles);
		
		int groupSize = ceil( (float)nTask / (float)nIOTasks );
		int groupMaster = round(Task.Rank / groupSize) * groupSize;
		int groupRank = Task.Rank - groupMaster;

		strncpy(filename, input_name, CHARBUFSIZE);

		if (restFiles >= nTask) { // read in nIO blocks, no communication
		
			int fileNum = nFiles - 1 - (Task.Rank + (restFiles - nTask)); 
			
			if (nFiles > 1)
				sprintf(filename, "%s.%i", filename, fileNum);

			for (int i = 0; i < groupSize; i++) {

				if (Task.Rank == groupMaster + i)
					read_file(filename, swapEndian, 0, 1, MPI_COMM_SELF);

				MPI_Barrier(MPI_COMM_WORLD);
			}

			restFiles -= nTask;

		} else { // parallel read on groupMaster  and distribute to group
		
			MPI_Comm_split(MPI_COMM_WORLD, groupMaster, groupRank, 
					&mpi_comm_read);
		
			int fileNum = nFiles - 1 - (restFiles - nIOTasks 
					+ groupMaster / groupSize);
			if (nFiles > 1)
				sprintf(filename, "%s.%i", filename, fileNum);

			read_file(filename, swapEndian, groupRank, groupSize, 
					mpi_comm_read);
			
			restFiles -= nIOTasks; 
			
			MPI_Comm_free(&mpi_comm_read);
		} 
		
		MPI_Barrier(MPI_COMM_WORLD);
	}
	
	if (RecvBuf != NULL)
		Free(RecvBuf); 
	if (ReadBuf != NULL)
		Free(ReadBuf);

	generate_masses_from_header();

	rprintf("\nReading completed\n\n");

	return ;
}

/* reads file on master and distributes it to an MPI communicator 
 * spanning groupSize with local rank groupRank */
static void read_file(char *filename, const bool swapEndian, 
		const int groupRank, const int groupSize, MPI_Comm mpi_comm_read)
{
	const int groupMaster = 0;  
	
	int i = 0, j = 0;

	FILE *fp = NULL;

	int nPartFile[NPARTYPE] = { 0 }, nTotRead = 0;

	if (groupRank == groupMaster) { // get nPart for this file

		fp = fopen(filename, "r");
	 	
		Assert(fp != NULL, "File not found %s", filename);

		find_block(fp, "HEAD", swapEndian);

		SKIP_FORTRAN_RECORD;
		
		safe_fread(nPartFile, sizeof(*nPartFile), 6, fp, swapEndian);
		
		SKIP_FORTRAN_RECORD;
		
		for (i=0; i<NPARTYPE; i++)
			nTotRead += nPartFile[i];
		
		printf("\nReading file '%s' on Task %i - %i \n"
        	    "   Gas   %9i   DM     %9i    \n"
           		"   Disk  %9i   Bulge  %9i    \n"
            	"   Star  %9i   Bndry  %9i    \n"
            	"   Total in File %9i \n\n",
				filename, Task.Rank, Task.Rank+groupSize-1,
				nPartFile[0], nPartFile[1], nPartFile[2],
				nPartFile[3], nPartFile[4], nPartFile[5], 
            	nTotRead); fflush(stdout);
	}

	MPI_Bcast(nPartFile, NPARTYPE, MPI_INT, 0, mpi_comm_read);
	
	int nPartGet[NPARTYPE] = { 0 }; // find npart per CPU and type
		
	for (i = 0; i < NPARTYPE; i++) 
		for (j = groupRank; j < nPartFile[i]; j += groupSize) 
			nPartGet[i]++;

	int nPartGetTotal = 0; 
	
	for (j = 0; j < NPARTYPE; j++) 
		nPartGetTotal += nPartGet[j];
	
	size_t nBytes = nPartGetTotal * largest_block_member_nbytes();
	RecvBuf = Realloc(RecvBuf, nBytes);

	nBytes = nTotRead * largest_block_member_nbytes(); // != 0 only for master
	ReadBuf = Realloc(ReadBuf, nBytes); 

	size_t offsets[NPARTYPE] = { 0 }; 

	Reallocate_P(nPartGet, offsets); // return offsets, update Task.Npart

	for (i = 0; i < NBlocks; i++) { // read blockwise

		int blocksize = 0;

		if (groupRank == groupMaster) { // find blocksize

			blocksize = find_block(fp, Block[i].Label, swapEndian);
		
			Assert(blocksize > 0 || !Block[i].PartBitMask, 
					"Can't find required block '%s'", Block[i].Label);

			printf("%18s %8d MB\n", Block[i].Name, blocksize/1024/1024);
			fflush(stdout);
		}

		MPI_Bcast(&blocksize, 1, MPI_INT, 0, mpi_comm_read);

		if (blocksize == 0) 
			continue ; // block not found
		
		if (groupRank == groupMaster) { // read on master
			
			nBytes = npart_in_block(i, nPartFile) * Block[i].Nbytes;
			
			Assert(nBytes == blocksize, 
					"File and Code blocksize inconsistent '%s', %zu != %d byte",
					Block[i].Label, nBytes, blocksize);

			SKIP_FORTRAN_RECORD
			
			safe_fread(ReadBuf, blocksize, 1, fp, swapEndian );
			
			SKIP_FORTRAN_RECORD
		} 
		
		int nBytesGetBlock = npart_in_block(i, nPartGet)* Block[i].Nbytes;
		
		int nBytesSend[groupSize];
		MPI_Allgather(&nBytesGetBlock, 1 , MPI_INT, nBytesSend, 1, MPI_INT, 
				mpi_comm_read);

		int displs[groupSize]; // displacements
		memset(displs, 0, groupSize*sizeof(*displs));

		for (j = 1; j < groupSize; j++) 
			displs[j] += nBytesSend[j-1] + displs[j-1];

		MPI_Scatterv(ReadBuf, nBytesSend, displs, MPI_BYTE, // distribute
					 RecvBuf, nBytesSend[groupRank], MPI_BYTE, 
					 groupMaster, mpi_comm_read); 

		empty_comm_buffer(RecvBuf, i, nPartGet, offsets); // copy *P 
	}

	if (fp != NULL)
		fclose(fp);
	
	return ;
}

static void empty_comm_buffer(char *DataBuf, const int iBlock, 
		const int *nPart, const size_t *offsets)
{
	const size_t nBytes = Block[iBlock].Nbytes;
	const size_t sizeof_P = sizeof(*P);

	char *start_P = (char *)P + Block[iBlock].Offset;
	//char *start_G = (char *)G + Block[iBlock].Offset;

	size_t src = 0; 

	switch (Block[iBlock].Target) { // ptr fun for the whole family

		case VAR_P:

			for (int type = 0; type < NPARTYPE; type++) { 
				
				size_t dest = offsets[type]*sizeof_P; 

				for (int i = 0; i < nPart[type]; i++) {  

					memcpy(start_P+dest, DataBuf+src, nBytes);
					
					dest += sizeof_P; // increments in bytes
					src += nBytes;
				}
			}
			
			break;

/*		case VAP_G:

			for (int i=0; i<nPart; i++)
				memcpy(&G[i]+offset, CommBuf+ nBytes*ibuf++, nBytes);

			break; */

		default: 
			Assert(0, "Input buffer target unknown");
	}

	return ;
}

static void read_header_data(FILE *fp, const bool swapEndian, int nFiles) 
{ 
	const bool swap = swapEndian;
	
	struct gadget_header head;

	int32_t blocksize = find_block(fp, "HEAD", swapEndian);

	Assert(blocksize == 256, "Format 2 Header corrupted");

	SKIP_FORTRAN_RECORD;

	safe_fread(head.Npart, sizeof(*head.Npart), 6, fp, swap);
	safe_fread(head.Massarr,sizeof(*head.Massarr), 6, fp, swap);
	safe_fread(&head.Time, sizeof(head.Time), 1, fp, swap);
	safe_fread(&head.Redshift, sizeof(head.Redshift), 1, fp, swap);
	safe_fread(&head.FlagSfr, sizeof(head.FlagSfr), 1, fp, swap);
	safe_fread(&head.FlagFeedback, sizeof(head.FlagFeedback), 1, fp,swap);
	safe_fread(head.Nall, sizeof(*head.Nall), 6, fp, swap);
	safe_fread(&head.FlagCooling, sizeof(head.FlagCooling), 1, fp,swap);
	safe_fread(&head.NumFiles, sizeof(head.NumFiles), 1, fp, swap);
	safe_fread(&head.Boxsize, sizeof(head.Boxsize), 1, fp, swap);
	safe_fread(&head.Omega0, sizeof(head.Omega0), 1, fp, swap);
	safe_fread(&head.OmegaLambda, sizeof(head.OmegaLambda), 1, fp, swap);
	safe_fread(&head.HubbleParam, sizeof(head.HubbleParam), 1, fp, swap);
	safe_fread(&head.FlagAge, sizeof(head.FlagAge), 1, fp, swap);
	safe_fread(&head.FlagMetals, sizeof(head.FlagMetals), 1, fp,swap);
	safe_fread(head.NallHighWord, sizeof(*head.NallHighWord), 6, fp,swap);

	for (int i= Sim.NpartTotal =0; i<NPARTYPE; i++) {
		Sim.Mpart[i] = head.Massarr[i];

		Sim.Npart[i] = (uint64_t)head.Nall[i];
		Sim.Npart[i] += ((uint64_t)head.NallHighWord[i]) << 32;

		Sim.NpartTotal += Sim.Npart[i];
	}

#ifdef PERIODIC
	Assert(Sim.Boxsize >= 0, "Boxsize in header not > 0, but %g ", Sim.Boxsize);
#endif	

	size_t sum = 0;

	for (int i=0; i<NPARTYPE; i++)
		sum += Sim.Mpart[i];

	printf("Total Particle Numbers (Masses) in Snapshot Header:	\n"
		"   Gas   %9llu (%1.5f), DM   %9llu (%1.5f), Disk %9llu (%1.5f)\n"
		"   Bulge %9llu (%1.5f), Star %9llu (%1.5f), Bndy %9llu (%1.5f)\n",
		(long long unsigned int)Sim.Npart[0], Sim.Mpart[0], 
		(long long unsigned int)Sim.Npart[1], Sim.Mpart[1], 
		(long long unsigned int)Sim.Npart[2], Sim.Mpart[2], 
		(long long unsigned int)Sim.Npart[3], Sim.Mpart[3], 
		(long long unsigned int)Sim.Npart[4], Sim.Mpart[4], 
		(long long unsigned int)Sim.Npart[5], Sim.Mpart[5]);
	
	if (head.NumFiles != nFiles)
		fprintf(stderr, "\nWARNING: NumFiles in Header (%d) "
				"doesnt match number of files for readin (%d) \n\n", 
				head.NumFiles, nFiles);

	return ;
}

static void generate_masses_from_header() 
{
	for (int i=0; i<NPARTYPE; i++)
		if (Sim.Mpart[i])
			return;

	int iMin = 0;
	
	for (int type=0; type<NPARTYPE; type++) {
		
		int iMax = iMin + Task.Npart[type];

		for (int ipart=iMin; ipart<iMax; ipart++)
			P[ipart].Mass = Sim.Mpart[type];

		iMin += Task.Npart[type];
	}

	return;
}

static int find_block(FILE *fp, const char label[4], const bool swapEndian)
{
	int32_t blocksize = 8;

	rewind(fp);

	for (;;) {

		char blocklabel[4] = {"    "};
		
		SKIP_FORTRAN_RECORD;
		
		safe_fread(blocklabel, 4 * sizeof(char), 1, fp, swapEndian);
		safe_fread(&blocksize, 4 * sizeof(char), 1, fp, swapEndian);

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

	return swapEndian;
}

static int find_files(char *filename)
{
	char buf[CHARBUFSIZE] = "";

	FILE *fp = fopen(filename, "r");

	int nFiles = 0;

	if (fp != NULL) { 

		nFiles = 1; 

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
static int safe_fread(void *data, size_t size, size_t nWanted, FILE *stream, 
		bool swapEndian)
{
	size_t nRead = fread(data, size, nWanted, stream);

	if (feof(stream)) // we catch EOF further up, because of block finding
		nRead = nWanted = 0;
		
	Assert(nRead == nWanted, "Read %zu bytes, but %zu wanted", nRead, nWanted);

	if (swapEndian && !feof(stream)) { // swap endianess

		char buf[16];

		for (int j=0; j<nWanted; j++) {

			memcpy(buf, &(( (char *)data )[j * size]), size);
			
			for (int i=0; i<size; i++)
				( (char*)data )[j*nWanted + i] = buf[size - i - 1];
		}
	}

	return nRead;
}

#undef SKIP_FORTRAN_RECORD
