#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>

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
