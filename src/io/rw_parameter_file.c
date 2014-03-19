#include "../globals.h"
#include "../proto.h"
#include "io.h"

void Read_Parameter_File(const char *filename)
{
	char buf[CHARBUFSIZE], buf1[CHARBUFSIZE],buf2[CHARBUFSIZE],
		 buf3[CHARBUFSIZE];
	bool tagDone[9999] = { false };
	
	if (!Task.Rank) {

		FILE *fd = fopen(filename, "r");

		Assert(fd != NULL, "Parameter file not found %s \n", filename);
			
		printf("\nReading Parameter file '%s' \n", filename);
			
		while (fgets(buf, CHARBUFSIZE, fd)) {

			sscanf(buf, "%s%s%s", buf1, buf2,buf3) ;

			if (buf1[0] == '%') // commented out
				continue;
			
			int j = -1;

			for (int i = 0; i < NTags; i++) {
				
				int tagNotFound = strcmp(buf1, ParDef[i].tag);
				
				if (!tagNotFound && !tagDone[i]) {

					j = i;

					tagDone[i] = true;
					
					break;
				}
			}

			if (j < 0) // don't know this one
				continue;
				
			printf(" %20s  %s\n", buf1, buf2);

			switch (ParDef[j].type) {

			case FLOAT:
				
				*((double *)ParDef[j].addr) = atof(buf2);
				
				break;
			case STRING:
			
				strncpy((char *)ParDef[j].addr, buf2, 
						CHARBUFSIZE);
				
				break;
			case INT:

				*((int *)ParDef[j].addr) = atoi(buf2);
				
				break;
			default:
				Assert(0, "Code Error in ParDef struct: %s",
						ParDef[j].tag);
			}
		}
		
		fclose(fd);
		
		printf("\n");

		for (int i = 0; i < NTags; i++) // are we missing one ?
            		Assert(tagDone[i],
				"Value for tag '%s' missing in parameter"
				" file '%s'.\n",ParDef[i].tag, filename );

	}

	MPI_Bcast(&Param, sizeof(Param), MPI_BYTE, 0, MPI_COMM_WORLD);

	/* sanity check parameters, yes yours too ! */
	Assert(Param.NumOutputFiles > 0, "NumOutputFiles has to be > 0");
	
	Assert(Param.NumIOTasks > 0, "NumIOTasks has to be > 0");
	
	Assert(Param.NumIOTasks <= Sim.NTask, 
		"NTask (=%d) can't be smaller than No_IOTasks (=%d)", 
		Sim.NTask,  Param.NumOutputFiles);
	
	Assert(Param.NumIOTasks <= Param.NumOutputFiles, 
		"NumIOTasks (=%d) can't be smaller than NumOutputFiles (=%d)", 
		Param.NumIOTasks,  Param.NumOutputFiles);
	
	Param.NumIOTasks = min(Param.NumIOTasks, Sim.NTask);
	
	Param.NumOutputFiles = min(Param.NumOutputFiles, Sim.NTask);
	
	MPI_Barrier(MPI_COMM_WORLD);

	return ;
}

void Write_Parameter_File(const char *filename)
{
	if (!Task.Rank) {

		printf("\nWriting Parameter file: %s \n", filename);

		FILE *fd = fopen(filename, "w");

		Assert(fd != NULL, "Can't open file %s for writing \n", filename);

		fprintf(fd, "%% Tandav, autogenerated parameter file %% \n\n");
		
		for (int i = 0; i < NTags; i++)
			fprintf(fd, "%s		%s \n", ParDef[i].tag, ParDef[i].val);

		fclose(fd);

		printf("\ndone, Good Bye.\n\n");
	}

	Finish();

	MPI_Finalize();

	exit(EXIT_SUCCESS);

	return ;
}
