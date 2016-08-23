/*!memory_management
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include <unistd.h>
#include <sys/stat.h>


int CreateDir(const char* pPath)
{
	if(-1 != access(pPath,0))
		return -1;

	char tmpPath[PATH_MAX];
	const char* pCur = pPath;

	memset(tmpPath,0,sizeof(tmpPath));
    
	int pos=0;
    while(*pCur++!='\0')
		{
			tmpPath[pos++] = *(pCur-1);

			if(*pCur=='/' || *pCur=='\0')
				{
					if(0!=access(tmpPath,0)&&strlen(tmpPath)>0)
						{
							mkdir(tmpPath, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
						}
				}
		}
	if(!access(pPath,0))
		return 0;
	else
		return 1;
}


void initialize_memory(double * p[],int N,int * CELL_POINT[])
{
	int k,i;
	for(k = 0; k < N; ++k)
		{
			p[k] = (double *)malloc(CELL_POINT[k][0] * sizeof(double));
			if(p[k] == NULL)
				{
					printf("Initialize memory fail.\n");
					exit(5);
				}
		}
}

void initialize_memory_int(int * p[],int N,int * CELL_POINT[])
{
	int k,i;
	for(k = 0; k < N; ++k)
		{
			p[k] = (int *)malloc(CELL_POINT[k][0] * sizeof(int));
			if(p[k] == NULL)
				{
					printf("Initialize memory fail.\n");
					exit(5);
				}
		}
}
