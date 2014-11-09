#include "../globals.h"
#include "../proto.h"
#include "../timestep.h"
#include "io.h"

#define WRITE_FORTRAN_RECORD(recSize) fwrite(&recSize, 4, 1, fp);

void write_file(const char *, const int, const int, const MPI_Comm);
void write_gadget_header(const int *npart, FILE *fp);
static void write_block_header(const char *, uint32_t, FILE *);
static void fill_data_buffer(const int, char *);
static void set_filename(char *filename);

static MPI_Comm mpi_comm_write = MPI_COMM_NULL;

void Write_Snapshot()
{ 
	Profile("Write Snap");

	#pragma omp single
	{

	const int nFiles = Param.Num_Output_Files;
	const int nIOTasks = Param.Num_IO_Tasks;

	int groupSize = Sim.NRank/nFiles; // big last file possible 
	int groupMaster = MIN(nFiles-1, floor(Task.Rank/groupSize)) * groupSize;

	int groupRank = Task.Rank - groupMaster;

	if (mpi_comm_write == MPI_COMM_NULL) // create & keep group communicator
		MPI_Comm_split(MPI_COMM_WORLD, groupMaster, groupRank, 
				&mpi_comm_write);

	int fileNum = groupMaster / groupSize;

	char filename[CHARBUFSIZE];

	set_filename(filename);

	if (nFiles > 1)
		sprintf(filename, "%s.%04i", filename, fileNum);
	
	for (int i = 0; i < nFiles; i+=nIOTasks) {

		if (fileNum < i+nIOTasks && fileNum >= i) 
			write_file(filename, groupRank, groupSize, mpi_comm_write);
		
		MPI_Barrier(MPI_COMM_WORLD);
	}

	Time.Snap_Counter++;

	} // omp single

	Profile("Write Snap");

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

	MPI_Reduce(&Task.Npart_Total, &nPartLargest, 1, MPI_INT, MPI_MAX, 
			groupMaster, mpi_comm_write);

	int nPartTotalFile = 0; // total number of particles in file

	for (int i = 0; i < NPARTYPE; i++) 
		nPartTotalFile += nPartFile[i];
	
	FILE *fp = NULL;
		
	if (groupRank == groupMaster) { // open file, write header

		printf("Writing file '%s' on MPI Ranks %i - %i \n"
			"   Gas   %9d, DM   %9d, Disk %9d\n"
			"   Bulge %9d, Star %9d, Bndy %9d\n"
			"   Total %9d\n\n",
				filename, Task.Rank, Task.Rank+groupSize-1,
				nPartFile[0],nPartFile[1],nPartFile[2],
				nPartFile[3],nPartFile[4],nPartFile[5],
				nPartTotalFile); 
		
		fp = fopen(filename, "w");
		
		Assert(fp != NULL, "Can't open file %s for writing", filename);
		
		write_gadget_header(nPartFile, fp);
	} 
	
	size_t dataBufSize = Largest_Block_Member_Nbytes(); 
	
	if (groupRank == groupMaster)
		dataBufSize *= 2*nPartLargest; // master stores comm&write buf 
	else
		dataBufSize *= Task.Npart_Total; // slaves buffer local data
	
	char *dataBuf = Malloc(dataBufSize, "dataBuf"); 
		
	for (int i = 0; i < NBlocks; i++) { // write blocks, hiding latency
	
		fill_data_buffer(i, dataBuf);

		size_t nBytesSend = Block[i].Nbytes * Npart_In_Block(i, Task.Npart);
		
		size_t xferSizes[groupSize]; // get size of data

		MPI_Gather(&nBytesSend, sizeof(nBytesSend), MPI_BYTE, 
					xferSizes,  sizeof(*xferSizes), 
					MPI_BYTE, groupMaster, mpi_comm_write);
	
		if (groupRank == groupMaster) { // master does all the work

			uint32_t blocksize = Npart_In_Block(i, nPartFile)
				* Block[i].Nbytes; 

			printf("%18s %8d MB\n", Block[i].Name, blocksize/1024/1024);

			write_block_header(Block[i].Label, blocksize, fp); 

			WRITE_FORTRAN_RECORD(blocksize)

			MPI_Request request; 
			MPI_Status status;
			
			int swap = 0; // to alternate between mem areas
			size_t halfBufSize = 0.5 * dataBufSize;

			/* recv & write */
			for (int task = 0; task < groupSize-1; task++) { 

				char * restrict writeBuf = dataBuf + swap * halfBufSize;
				char * restrict commBuf = dataBuf + (1 - swap) * halfBufSize;
		
				MPI_Irecv(commBuf, xferSizes[task+1], MPI_BYTE,	task+1, task+1,
						mpi_comm_write, &request); 

				fwrite(writeBuf, xferSizes[task], 1, fp);

				swap = 1 - swap; // swap memory areas 

				MPI_Wait(&request, &status);
			}

			/* last one in group */
			fwrite(dataBuf+swap*halfBufSize, xferSizes[groupSize-1], 1, fp);

			WRITE_FORTRAN_RECORD(blocksize) 

		} else   // slaves just post a blocking send
			MPI_Send(dataBuf, nBytesSend, MPI_BYTE, groupMaster, 
					groupRank, mpi_comm_write);

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
	
	uint32_t blocksize = sizeof(head);
	
	Assert(blocksize == 256, "sizeof(head) incorrect, %d byte", blocksize);

	for (int i = 0; i < 6; i++) {
		
		head.Npart[i] = npart[i];

		head.Nall[i] = (int32_t)(Sim.Npart[i]);
		head.Nall_High_Word[i] = (int32_t) (Sim.Npart[i] >> 32);

		head.Massarr[i] = Sim.Mpart[i];
	}

	head.Time = Time.Current;

#ifdef COMOVING
	head.Redshift = 1.0/(1.0 + head.Time);
#endif // COMOVING

	head.Flag_Sfr = 0;
	head.Flag_Feedback = 0;
	head.Flag_Cooling = 0;
	head.Num_Files = Param.Num_Output_Files;
	head.Boxsize = Sim.Boxsize[0]; // fall back 
	head.Omega0 = Cosmo.Omega_0;
	head.Omega_Lambda = Cosmo.Omega_Lambda;
	head.Hubble_Param = Cosmo.Hubble_Constant;
	head.Flag_Age = 0;
	head.Flag_Metals = 0;

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

			src = (char *) P + offset;
			
			for (int i = 0; i < Task.Npart_Total; i++) {
			
				memcpy(dest, src, nBytes);

				dest += nBytes;
				src += sizeof_P;
			}
			
			break;

		case VAR_GAS:
			break;

		case VAR_DM:
			break;

		case VAR_STAR:
			break;

		case VAR_DISK:
			break;

		case VAR_BND:
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

static void set_filename(char *filename)
{
	const int ndigits = ceil(log10(Time.NSnap));

	switch (ndigits) {
		
	case 4:
		sprintf(filename, "%s_%04d", 
				Param.Output_File_Base, Time.Snap_Counter);
		break;

	case 5:
		sprintf(filename, "%s_%05d", 
			Param.Output_File_Base, Time.Snap_Counter);
		break;
		
	default:
		sprintf(filename, "%s_%03d", 
				Param.Output_File_Base, Time.Snap_Counter);
		break;
	}

	return;
}
