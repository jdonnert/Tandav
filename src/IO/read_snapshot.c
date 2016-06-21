#include "io.h"

#define SKIP_FORTRAN_RECORD safe_fread(&Fortran_Record, 4, 1, fp, swap_Endian);

static int32_t Fortran_Record; // holds the 4 byte Fortran data

static int safe_fread(void * restrict, size_t, size_t, FILE *, bool);
static int find_files(char *);
static bool find_endianess(FILE *);
static int find_block(FILE *, const char label[4], const bool);
static void read_header_data(FILE *fp, const bool, int);
static void read_file(char *, const bool, const int, const int, MPI_Comm);
static void empty_comm_buffer(char * restrict, const int, const int, 
		const int *, const size_t *);

static void generate_masses_from_header();
static void set_particle_types();

void Read_Snapshot(char *input_name)
{
	int nIOTasks = Param.Num_IO_Tasks;

	int nFiles = 0;
	bool swap_Endian = false;
	char filename[CHARBUFSIZE] = "";

	Profile("Read Snap");

	if (Task.Is_Master) {

		nFiles = find_files(input_name);

		rprintf("\nParallel Reading of %d files on %d tasks\n\n",
				nFiles, nIOTasks);

		strncpy(filename, input_name, CHARBUFSIZE);

		if (nFiles > 1)
			sprintf(filename, "%s.0", input_name);

		FILE *fp = fopen(filename, "r");

		swap_Endian = find_endianess(fp);

		read_header_data(fp, swap_Endian, nFiles); // fills Sim

		fclose(fp);
	}
	
	MPI_Bcast(&nFiles, 1, MPI_INT, MASTER, MPI_COMM_WORLD);

	MPI_Bcast(&swap_Endian, sizeof(swap_Endian), MPI_BYTE, MASTER, 
			MPI_COMM_WORLD);

	MPI_Bcast(&Sim, sizeof(Sim), MPI_BYTE, MASTER, MPI_COMM_WORLD);

	Allocate_Particle_Structures();
	
	int rest_Files = nFiles;

	while (rest_Files > 0) {

		nIOTasks = fmin(nIOTasks,rest_Files);
		
		int groupSize = ceil( (double)NRank / (double)nIOTasks );
		int groupMaster = round(Task.Rank / groupSize) * groupSize;
		int groupRank = Task.Rank - groupMaster;

		strncpy(filename, input_name, CHARBUFSIZE);

		if (rest_Files >= NRank) { // read in nIO blocks, no communication
		
			int fileNum = nFiles - 1 - (Task.Rank + (rest_Files - NRank)); 
			
			if (nFiles > 1)
				sprintf(filename, "%s.%i", filename, fileNum);

			for (int i = 0; i < groupSize; i++) {

				if (Task.Rank == groupMaster + i)
					read_file(filename, swap_Endian, 0, 1, MPI_COMM_SELF);

				MPI_Barrier(MPI_COMM_WORLD);
			}

			rest_Files -= NRank;

		} else { // parallel read on groupMasters and distribute to group
	
			MPI_Comm mpi_comm_read;
		
			MPI_Comm_split(MPI_COMM_WORLD, groupMaster, groupRank, 
					&mpi_comm_read);
		
			int fileNum = nFiles - 1 - (rest_Files - nIOTasks 
					+ groupMaster / groupSize);

			if (nFiles > 1)
				sprintf(filename, "%s.%i", filename, fileNum);

			read_file(filename, swap_Endian, groupRank, groupSize, 
					mpi_comm_read);
			
			rest_Files -= nIOTasks; 
			
			MPI_Comm_free(&mpi_comm_read);
		} 
		
		MPI_Barrier(MPI_COMM_WORLD);
	}

	generate_masses_from_header();

 	set_particle_types();

	rprintf("\nReading completed\n\n");

	#pragma omp barrier

	Profile("Read Snap");

	return ;
}

/* 
 * reads file on master and distributes it to an MPI communicator 
 * spanning groupSize with local rank groupRank 
 */

static void read_file (char *filename, const bool swap_Endian, 
		const int groupRank, const int groupSize, MPI_Comm mpi_comm_read)
{
	const int groupMaster = 0;  
	
	int i = 0, j = 0;

	FILE *fp = NULL;

	int nPartFile[NPARTYPE] = { 0 }, nTotRead = 0;

	if (groupRank == groupMaster) { // get nPart for this file

		fp = fopen(filename, "r");
	 	
		Assert(fp != NULL, "File not found %s", filename);

		find_block(fp, "HEAD", swap_Endian);

		SKIP_FORTRAN_RECORD;
		
		safe_fread(nPartFile, sizeof(*nPartFile), 6, fp, swap_Endian);
		
		SKIP_FORTRAN_RECORD;
		
		for (i = 0; i < NPARTYPE; i++)
			nTotRead += nPartFile[i];
		
		printf("\nReading file '%s' on Task %i - %i \n"
       	   	"   Gas   %9i   DM     %9i    \n"
       		"   Disk  %9i   Bulge  %9i    \n"
           	"   Star  %9i   Bndry  %9i    \n"
           	"   Total in File %9i \n\n",
			filename, Task.Rank, Task.Rank+groupSize-1,
			nPartFile[0], nPartFile[1], nPartFile[2],
			nPartFile[3], nPartFile[4], nPartFile[5], 
       		nTotRead); 
		fflush(stdout);
	}

	MPI_Bcast(nPartFile, NPARTYPE, MPI_INT, 0, mpi_comm_read);
	
	int nPartGet[NPARTYPE] = { 0 }; // find npart per CPU and type
		
	for (i = 0; i < NPARTYPE; i++) 
		for (j = groupRank; j < nPartFile[i]; j += groupSize) 
			nPartGet[i]++;
	
	size_t offsets[NPARTYPE] = { 0 }; 
	
	Reallocate_P(nPartGet, offsets); // return offsets, update Task.Npart

	int nPartGetTotal = 0; 
	
	for (j = 0; j < NPARTYPE; j++) 
		nPartGetTotal += nPartGet[j];
	
	size_t nBytes = nPartGetTotal * Largest_Block_Member_Nbytes();

	char *RecvBuf = Malloc(nBytes, "RecvBuf");

	char *ReadBuf = NULL;

	if (groupRank == groupMaster) {
	
		nBytes = nTotRead * Largest_Block_Member_Nbytes(); 
		
		ReadBuf = Malloc(nBytes, "ReadBuf"); 
	}

	for (i = 0; i < NBlocks; i++) { // read blockwise

		uint32_t blocksize = 0;

		if (groupRank == groupMaster) { // find blocksize

			blocksize = find_block(fp, Block[i].Label, swap_Endian);
		
			Assert(blocksize != 0 || (Block[i].IC_Required == false), 
					"Can't find required block '%s'", Block[i].Label);
			
			if (blocksize > 0)
				printf("%18s %8d MB\n", Block[i].Name, blocksize/1024/1024);

			fflush(stdout);
		}
		
		MPI_Bcast(&blocksize, 1, MPI_INT, 0, mpi_comm_read);

		if (blocksize == 0) 
			continue ; // block not found
		
		if (groupRank == groupMaster) { // read on master
			
			nBytes = Npart_In_Block(i, nPartFile) * Block[i].Ncomp 
					* Block[i].Nbytes;
			
			Assert(nBytes == blocksize, 
				"File and Code blocksize inconsistent '%s', %zu != %d byte", 
				Block[i].Label, nBytes, blocksize);

			SKIP_FORTRAN_RECORD
			
			safe_fread(ReadBuf, nBytes, 1, fp, swap_Endian);

			SKIP_FORTRAN_RECORD
		
		} 
		
		int nBytesGetBlock = Npart_In_Block(i, nPartGet) 
										* Block[i].Ncomp * Block[i].Nbytes;
		
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

		empty_comm_buffer(RecvBuf, i, nPartGetTotal, nPartGet, offsets);  
	}

	if (fp != NULL)
		fclose(fp);
	
	Free(RecvBuf); Free(ReadBuf);
	
	return ;
}

/* 
 * This moves data from the comm buffer to P. Here we should be limited 
 * by the I/O of the drive. We find the destination of the block by counting
 * the number of pointers (nPtr) starting from &P.Type. nPtr is 
 * just the offset in Bytes over the size of a pointer in bytes on this 
 * system.
 */

static void empty_comm_buffer(char * restrict DataBuf, const int iB, 
		const int nPartTotal, const int *nPart, const size_t *offsets)
{
	const int nComp = Block[iB].Ncomp;
	const size_t nBytes = Block[iB].Nbytes;
	const size_t nPtr = Block[iB].Offset/sizeof(void *); 
			
	char * restrict src = DataBuf; 
	char * restrict dest[nComp];

	switch (Block[iB].Target) { // find the destination pointers

		case VAR_P:

			for (int j = 0; j < nComp; j++) // ptr fun for the whole family
				dest[j] = (char *) *(&P.Type + nPtr + j); // points to P.X[i]

		break;

		default: 

			Assert(false, "Input buffer target unknown %d", Block[iB].Target);

		break;
	}	
	
	for (int i = 0; i < nPartTotal; i++) {
			
		for (int j = 0; j < nComp; j++) {
				
			memcpy(dest[j], src, nBytes);
				
			src += nBytes;
			dest[j] += nBytes;
		}
	}

	return ;
}

static void read_header_data(FILE *fp, const bool swap_Endian, int nFiles) 
{ 
	const bool swap = swap_Endian;
	
	struct gadget_header head = { { 0 } };

	int32_t blocksize = find_block(fp, "HEAD", swap_Endian);

	Assert(blocksize == 256, "Format 2 Header corrupted");

	SKIP_FORTRAN_RECORD;

	safe_fread(head.Npart, sizeof(*head.Npart), 6, fp, swap);
	safe_fread(head.Massarr,sizeof(*head.Massarr), 6, fp, swap);
	safe_fread(&head.Time, sizeof(head.Time), 1, fp, swap);
	safe_fread(&head.Redshift, sizeof(head.Redshift), 1, fp, swap);
	safe_fread(&head.Flag_Sfr, sizeof(head.Flag_Sfr), 1, fp, swap);
	safe_fread(&head.Flag_Feedback, sizeof(head.Flag_Feedback), 1, fp,swap);
	safe_fread(head.Nall, sizeof(*head.Nall), 6, fp, swap);
	safe_fread(&head.Flag_Cooling, sizeof(head.Flag_Cooling), 1, fp,swap);
	safe_fread(&head.Num_Files, sizeof(head.Num_Files), 1, fp, swap);
	safe_fread(&head.Boxsize, sizeof(head.Boxsize), 1, fp, swap);
	safe_fread(&head.Omega0, sizeof(head.Omega0), 1, fp, swap);
	safe_fread(&head.Omega_Lambda, sizeof(head.Omega_Lambda), 1, fp, swap);
	safe_fread(&head.Hubble_Param, sizeof(head.Hubble_Param), 1, fp, swap);
	safe_fread(&head.Flag_Age, sizeof(head.Flag_Age), 1, fp, swap);
	safe_fread(&head.Flag_Metals, sizeof(head.Flag_Metals), 1, fp,swap);
	safe_fread(head.Nall_High_Word, sizeof(*head.Nall_High_Word), 6, fp,swap);

	Sim.Npart_Total = 0;

	for (int i = 0; i < NPARTYPE; i++) {

		Sim.Mpart[i] = head.Massarr[i];

		Sim.Npart[i] = (uint64_t)head.Nall[i];
		Sim.Npart[i] += ((uint64_t)head.Nall_High_Word[i]) << 32;
		
		Sim.Npart_Total += Sim.Npart[i];
	}

	Assert(head.Boxsize >= 0, "Boxsize in header not > 0, but %g", Sim.Boxsize);

#ifdef PERIODIC
	if (Sim.Boxsize[0] == -1) {

		Sim.Boxsize[0] = Sim.Boxsize[1] = Sim.Boxsize[2] = head.Boxsize;

		printf("Setting boxsize from snapshot header: %g \n\n",Sim.Boxsize[0]);
	}

#endif	

	size_t sum = 0;

	for (int i = 0; i < NPARTYPE; i++)
		sum += Sim.Npart[i];

	printf("Snapshot Header:"
		" %d Files, Boxsize of %g at time %g \n\n"
		"Total Particle Numbers (Masses):	\n"
		"   Gas   %9llu (%5.2g), DM   %9llu (%5.2g), Disk %9llu (%5.2g)\n"
		"   Bulge %9llu (%5.2g), Star %9llu (%5.2g), Bndy %9llu (%5.2g)\n"
		"   Sum %10zu \n",
		head.Num_Files, head.Boxsize, head.Time,
		(long long unsigned int) Sim.Npart[0], Sim.Mpart[0], 
		(long long unsigned int) Sim.Npart[1], Sim.Mpart[1], 
		(long long unsigned int) Sim.Npart[2], Sim.Mpart[2], 
		(long long unsigned int) Sim.Npart[3], Sim.Mpart[3], 
		(long long unsigned int) Sim.Npart[4], Sim.Mpart[4], 
		(long long unsigned int) Sim.Npart[5], Sim.Mpart[5], sum);

	Assert(head.Num_Files == nFiles, "NumFiles in Header (%d) doesnt match "
			"number of files found (%d) \n\n", head.Num_Files, nFiles);

	Warn(head.Omega0 != Cosmo.Omega_Matter,
			"Omega_0 in snapshot different from code: file %g, code %g",
			head.Omega0, Cosmo.Omega_Matter);

	Warn(head.Omega_Lambda != Cosmo.Omega_Lambda,
			"Omega_Lambda in snapshot different from code: file %g, code %g",
			head.Omega_Lambda, Cosmo.Omega_Lambda);

	Warn(head.Hubble_Param != HUBBLE_CONST/100.0,
			"h_0 in snapshot different from code: file %g, code  %g",
			head.Hubble_Param, HUBBLE_CONST/100.0);

	Warn(head.Boxsize != Sim.Boxsize[0], "Boxsize inconsistent %g <-> %g,%g,%g",
			head.Boxsize, Sim.Boxsize[0], Sim.Boxsize[1], Sim.Boxsize[2]);

	if (Param.Start_Flag == READ_SNAP) 
		Restart.Time_Continue = head.Time;

	return ;
}

static void generate_masses_from_header()
{
	bool mass_in_header = false;

	for (int type = 0; type < NPARTYPE; type++)
		if (Sim.Mpart[type] != 0)
			mass_in_header = true;

	if (! mass_in_header)
		return;

	int iMin = 0;

	for (int type = 0; type < NPARTYPE; type++) {

		int iMax = iMin + Task.Npart[type];
		
		#pragma  omp parallel for shared(iMin,iMax)
		for (int ipart = iMin; ipart < iMax; ipart++)
			P.Mass[ipart] = Sim.Mpart[type];

		iMin += Task.Npart[type];
	}

	return;
}

/*
 * We recover the particle types from the IDs which are assumed strictly 
 * ordered in the range [1, Npart]
 */

static void set_particle_types()
{
	uint64_t min_ID[NPARTYPE] = { 0 };
	uint64_t max_ID[NPARTYPE] = { 0 };

	uint64_t run = 0;
		
	for (int type = 0; type < NPARTYPE; type++) {

		min_ID[type] = max_ID[type] = -1;

		if (Sim.Npart[type] == 0)
			continue;

		min_ID[type] = run;
		
		run += Sim.Npart[type];

		max_ID[type] = run - 1;
	}

	#pragma omp parallel for shared(max_ID,min_ID)
	for (int ipart = 0; ipart < Task.Npart_Total; ipart++) {
		
		for (int type = 0; type < NPARTYPE; type++) {
			
			if (P.ID[ipart] <= max_ID[type])
				if (P.ID[ipart] >= min_ID[type])
					P.Type[ipart] = type;
		}
	}
	
	return;
}

static int find_block(FILE *fp, const char label[4], const bool swap_Endian)
{
	int32_t blocksize = 8;

	rewind(fp);

	for (;;) {

		char blocklabel[4] = {"    "};
		
		SKIP_FORTRAN_RECORD;
		
		safe_fread(blocklabel, 4 * sizeof(char), 1, fp, swap_Endian);
		safe_fread(&blocksize, 4 * sizeof(char), 1, fp, swap_Endian);

        SKIP_FORTRAN_RECORD;

        if (strncmp(label, blocklabel, 4) == 0) 
			break; // found it
		
		if (!feof(fp)) {
		
			fseek(fp, blocksize, SEEK_CUR); // skip to end of block
		
		} else {

			blocksize = 8;

			break;
		}
	}

	return blocksize - 8; // remove 8 bytes from Fortran header
}


static bool find_endianess(FILE *fp) 
{
	rewind(fp);

	int32_t testRecord = 0;
	
	size_t nRead = fread(&testRecord, 4, 1, fp);

	Assert(nRead == 1, "Couldn't read Fortran test record for Endianess");

	bool swap_Endian = false;
	
	if (testRecord == 134217728) { // first block is always 8 bytes

		printf("\nEnabling Endian Swapping\n");

		swap_Endian = true;
	}

	rewind(fp);

	SKIP_FORTRAN_RECORD

	Assert(Fortran_Record == 8, "Binary Fortran File Format Broken");

	return swap_Endian;
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

            Assert(nFiles < 10000, "Found 10000 files, holy cow !");
		}

    	Assert(nFiles, "Can't open input file as '%s' or '%s'", filename, buf);
	}

	return nFiles;
}

/* 
 * this fread wrapper checks for eof, corruption and swaps endianess 
 * if swap_Endian == true 
 */

static int safe_fread(void * restrict data, size_t size, size_t nWanted, 
		FILE *stream, bool swap_Endian)
{
	size_t nRead = fread(data, size, nWanted, stream);

	if (feof(stream)) // we catch EOF further up, because of block finding
		nRead = nWanted = 0;
	
	Assert(nRead == nWanted, "Read %zu blocks with %zu bytes, but %zu wanted",
			nRead, size, nWanted);

	if (swap_Endian && !feof(stream)) { // swap endianess

		char buf[16] = { " " };

		for (int j = 0; j < nWanted; j++) {

			memcpy(buf, &(( (char *)data )[j * size]), size);
			
			for (int i = 0; i < size; i++) {

				size_t dest = j*nWanted + i;
				
				size_t src = size - i - 1;

				((char*) data)[dest] = buf[src];
			}
		}
	}

	return nRead;
}

#undef SKIP_FORTRAN_RECORD
