#include "../globals.h"
#include "io.h"

static size_t save_fread(void *, size_t size, size_t nWanted, FILE *stream);
static void read_file(char *filename, int ReadTask, int LastTask);
static int find_files(char *filename);
static bool find_endianess(FILE * fp);
static void read_header_data(FILE *fp);

#define READ_FORTRAN_RECORD save_fread(&FortranRecord, 4, 1, fp, SwitchEndian);
static int32_t FortranRecord; // holds the 4 byte Fortran header block

static bool SwitchEndian = false;
static int NFiles = 0;

void Read_Snapshot(char *filename)
{
	if (!Task.Rank) {
		NFiles = find_files(filename);
		
		if (NFiles > 1)
    	    sprintf(filename, "%s.0", filename);

		FILE *fp = fopen(filename, "r");
		Assert(fp != NULL, "File %s disappeared !",filename);
		
		SwitchEndian = find_endianess(fp);

		read_header_data(fp); // fills global Sim variable
	}

	MPI_Bcast(&NFiles, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&SwitchEndian, sizeof(SwitchEndian), MPI_BYTE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&Sim, sizeof(Sim), MPI_BYTE, 0, MPI_COMM_WORLD);

	const int nIOTasks = Param.No_Output_Files;
	rprintf("\nParallel Reading on %d tasks\n",	nIOTasks);




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

static void read_header_data(FILE *fp) { 

	return ;
}

static bool find_endianess(FILE *fp) 
{
	rewind(fp);

	int32_t testRecord = 0;
	
	size_t nRead = fread(&testRecord, 4, 1, fp);
	Assert(nRead == 1, "Couldn't read Fortran record");

	if (testRecord == 134217728) { // first block is always 8 bytes in format 2
		SwitchEndian = true;
		printf("\nEnabling Endian Swapping\n");
	}

	rewind(fp);

	save_fread(&FortranRecord, 4, 1, fp); // if needed this now switches endian

	Assert(FortranRecord == 8, "First Fortran record != 8");

	fclose(fp);

	return SwitchEndian;
}

int find_files(char *filename)
{
	char buf[MAXSTRINGLENGTH] = "\0";

	FILE *fp = fopen(filename, "r");

	if (fp != NULL) { 
		NFiles++;

		fclose(fp);
	} else {
		for (;;) {
            sprintf(buf, "%s.%i", filename, NFiles);
			
            if (!(fp = fopen(buf, "r"))) 
				break;
			
			fclose(fp);
			
            NFiles++;

            Assert(NFiles<1000, "Found more than 1000 files ??");
		}

    	Assert(NFiles, "Can't open input file as '%s' or '%s'", 
				filename, buf);
	}

	return NFiles;
}

static size_t save_fread(void *data, size_t size, size_t nWanted, FILE *stream)
{
	size_t nRead = fread(data, size, nWanted, stream);

	if (feof(stream)) // we catch EOF further up, because of block finding
		nRead = nWanted = 0;
		
	Assert(nRead == nWanted, "Read %zu bytes, but %zu wanted", nRead, nWanted);

	if (SwitchEndian) { // swap endianess
		char buf[16];

		for (int j=0; j < nWanted; j++) {
			memcpy(buf, &( ( (char *)data) [j * size]), size);
			
			for (int i=0; i < size; i++)
				( (char*)data )[j*nWanted + i] = buf[size - i - 1];
		}
	}

	return nRead;
}
