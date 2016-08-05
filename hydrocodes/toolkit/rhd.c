#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double rnd( double *r)//random number
{
    int m;
    double s,u,v,p;
    s=65536.0; u=2053.0; v=13849.0;
    m=(int)(*r/s); *r=*r-m*s;
    *r=u*(*r)+v; m=(int)(*r/s);
    *r=*r-m*s; p=*r/s;
    return(p);
}

