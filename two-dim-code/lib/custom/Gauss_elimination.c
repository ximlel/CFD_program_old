#include <stdio.h>
#include <math.h>
#include <stdlib.h>


void Gauss_elimination(int n, double (*a)[n+1], double *x)
{ 

	int i,j,k;
	double temp,s,l;

	for(i=0;i<n-1;i++)
		{
			//selected principal element of columns
			k=i;     
			for(j=i+1;j<n;j++)
				{
					if(fabs(a[j][i])>fabs(a[k][i]))
						k=j;
				}
			
			//line feed   
			if(k!=i)
				for(j=i;j<=n;j++)
					{
						temp=a[i][j];
						a[i][j]=a[k][j];
						a[k][j]=temp;
					}
			
			//elimination
			for(j=i+1;j<n;j++)
				{ 
					l=a[j][i]/a[i][i];
					for(k=0;k<n+1;k++)
						a[j][k]=a[j][k]-a[i][k]*l;      
				}			
		}
	//back substitution
	x[n-1]=a[n-1][n]/a[n-1][n-1];
	for(i=n-2;i>=0;i--)
		{
			s=0.0;
			for(j=i;j<n;j++)
				{
					if(j==i)
						continue;
					s+=a[i][j]*x[j];					
				}
			x[i]=(a[i][n]-s)/a[i][i];
		}
}
