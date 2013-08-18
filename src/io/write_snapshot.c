#include "../globals.h"
#include "io.h"

#define WRITE_FORTRAN_RECORD(recSize) fwrite(&recSize, 4, 1, fp);

void write_file(const char *, const int, const int, const MPI_Comm);
void write_gadget_header(const int *npart, FILE *fp);
static void write_block_header(const char *, uint32_t, FILE *);
static void fill_data_buffer(const int, char *);

void Write_Snapshot()
{
	const int nFiles = Param.NumOutputFiles;
	const int nIOTasks = Param.NumIOTasks;
	
	MPI_Comm mpi_comm_write;

	int groupSize = Sim.NTask/nFiles; // big last file possible 
	int groupMaster = min(nFiles-1, floor(Task.Rank/groupSize)) * groupSize;

	int groupRank = Task.Rank - groupMaster;

	MPI_Comm_split(MPI_COMM_WORLD, groupMaster, groupRank, &mpi_comm_write);

	int fileNum = groupMaster / groupSize;

	char filename[CHARBUFSIZE];

	sprintf(filename, "%s_%03d",Param.OutputFileBase, Time.SnapCounter);
	
	if (nFiles > 1)
		sprintf(filename, "%s.%i", filename, fileNum);
	
	for (int i = 0; i < nFiles; i+=nIOTasks) {

		if (fileNum < i+nIOTasks && fileNum >= i) 
			write_file(filename, groupRank, groupSize, mpi_comm_write);
		
		MPI_Barrier(MPI_COMM_WORLD);
	}

	Time.SnapCounter++;

	rprintf("Snapshot written\n");

	return ;
}

void write_file(const char *filename, const int groupRank, const int groupSize,
		const MPI_Comm mpi_comm_write)
{
	const int groupMaster = 0;  

	int nPartFile[NPARTYPE] = { 0 }; // npart in file by type
	MPI_Reduce(Task.Npart, nPartFile, NPARTYPE, MPI_INT, MPI_SUM, 
			groupMaster, mpi_comm_write);

	int nPartLargest = 0; // largest part number to transfer
	MPI_Reduce(&Task.NpartTotal, &nPartLargest, 1, MPI_INT, MPI_MAX, 
			groupMaster, mpi_comm_write);

	int nPartTotalFile = 0; // total number of particles in file
	for (int i = 0; i < NPARTYPE; i++) 
			nPartTotalFile += nPartFile[i];
	
	FILE *fp = NULL;
		
	if (groupRank == groupMaster) { // open file, write header

		printf("Writing file '%s' on Task %i - %i \n"
			"   Gas   %9d, DM   %9d, Disk %9d\n"
			"   Bulge %9d, Star %9d, Bndy %9d\n"
			"   Total %9d\n\n",
				filename, Task.Rank, Task.Rank+groupSize-1,
				nPartFile[0],nPartFile[1],nPartFile[2],
				nPartFile[3],nPartFile[4],nPartFile[5],nPartTotalFile);

		fp = fopen(filename, "w");
		
		Assert(fp != NULL, "Can't open file %s for writing", filename);
		
		write_gadget_header(nPartFile, fp);
	} 
	
	size_t dataBufSize = largest_block_member_nbytes(); 
	
	if (groupRank == groupMaster)
		dataBufSize *= 2*nPartLargest; // master stores comm & write buf 
	else
		dataBufSize *= Task.NpartTotal; // slaves just local data
	
	char *dataBuf = Malloc(dataBufSize); 
		
	for (int i = 0; i < NBlocks; i++) { // write blocks, hiding comm. latency

		fill_data_buffer(i, dataBuf);

		size_t nBytesSend = Block[i].Nbytes * npart_in_block(i, Task.Npart);
		
		size_t xferSizes[groupSize]; // get size of data

		MPI_Gather(&nBytesSend, sizeof(nBytesSend), MPI_BYTE, 
					xferSizes,  sizeof(*xferSizes), MPI_BYTE, 
					groupMaster, mpi_comm_write);
	
		if (groupRank == groupMaster) { // master does all the work

			uint32_t blocksize = npart_in_block(i, nPartFile)*Block[i].Nbytes; 

			write_block_header(Block[i].Label, blocksize, fp); 
			
			WRITE_FORTRAN_RECORD(blocksize)

			MPI_Request request; MPI_Status status;
			
			size_t swap = 0; // to swap memory areas
			size_t halfBufSize = 0.5 * dataBufSize;

			for (int task = 0; task < groupSize-1; task++) {

				char * restrict writeBuf = dataBuf + swap * halfBufSize;
				char * restrict commBuf = dataBuf + (1 - swap) * halfBufSize;
		
				MPI_Irecv(commBuf, xferSizes[task+1], MPI_BYTE, task+1, 
						task+1, mpi_comm_write, &request); 

				fwrite(writeBuf, xferSizes[task], 1, fp);

				swap = 1 - swap; // swap memory areas 

				MPI_Wait(&request, &status);
			}

			/* last one in group */
			fwrite(dataBuf+swap*halfBufSize, xferSizes[groupSize-1], 1, fp);

			WRITE_FORTRAN_RECORD(blocksize) 

		} else   // slaves just post a blocking send
			MPI_Send(dataBuf, nBytesSend, MPI_BYTE, groupMaster, groupRank, 
					mpi_comm_write);

		MPI_Barrier(mpi_comm_write);
	}

	if (groupRank == groupMaster)
		fclose(fp);
	
	Free(dataBuf);

	MPI_Barrier(mpi_comm_write);

	return ;
}

void write_gadget_header(const int *npart, FILE *fp) 
{
	struct gadget_header head;
	
	Assert(sizeof(head)==256, "sizeof(head) incorrect, %d byte", sizeof(head));

	for (int i = 0; i < 6; i++) {
		
		head.Npart[i] = npart[i];

		head.Nall[i] = (int32_t)(Sim.Npart[i]);
		head.NallHighWord[i] = (int32_t) (Sim.Npart[i] >> 32);

		head.Massarr[i] = Sim.Mpart[i];
	}

	head.Time = Time.Current;
	head.Redshift = 1.0/(1.0+head.Time);
	head.FlagSfr = 0;
	head.FlagFeedback = 0;
	head.FlagCooling = 0;
	head.NumFiles = Param.NumOutputFiles;
	head.Boxsize = Sim.Boxsize;
	head.Omega0 = Cosmo.Omega0;
	head.OmegaLambda = Cosmo.OmegaLambda;
	head.HubbleParam = Cosmo.HubbleParam;
	head.FlagAge = 0;
	head.FlagMetals = 0;

	uint32_t blocksize = sizeof(head);
	
	write_block_header("HEAD", blocksize, fp);
		
	WRITE_FORTRAN_RECORD(blocksize)

	fwrite(&head, blocksize, 1, fp);
	
	WRITE_FORTRAN_RECORD(blocksize)

	return ;
}

static void fill_data_buffer(const int i, char *dataBuf)
{
	const size_t sizeof_P = sizeof(*P);
	//const size_t sizeof_G = sizeof(*G);
	
	const size_t offset = Block[i].Offset;
	const size_t nBytes = Block[i].Nbytes;

	char * restrict src = NULL;
	char * restrict dest = dataBuf;

	switch (Block[i].Target) {

		case VAR_P:

			src = (char *) P+offset;
			
			for (int i = 0; i < Task.NpartTotal; i++) 
				memcpy(dest+i*nBytes, src+i*sizeof_P, nBytes);
			
			break;

		case VAR_G:
			break;

		default:
			Assert(0, "Block Target not handled : %d", Block[i].Target);
	}

	return ;
}
static void write_block_header(const char *name, uint32_t blocksize, FILE *fp) 
{
	Assert(blocksize <= UINT_MAX - 8, 
			"Block %s too large to fit FORTRAN format", name);

	const uint32_t fmt2Size = 8; // size of the format 2 header

	blocksize += 8; // add 2*4 byte of FORTRAN header to data size

	char fmt2Head[fmt2Size];

	strncpy(&fmt2Head[0], name, 4);
	strncpy(&fmt2Head[4], (char *)&blocksize, 4);

	WRITE_FORTRAN_RECORD(fmt2Size)

	fwrite(&fmt2Head, fmt2Size, 1, fp);
	
	WRITE_FORTRAN_RECORD(fmt2Size)

	return ;
}
