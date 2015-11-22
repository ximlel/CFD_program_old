double miu_BJ(double x) //Barth_Jespersen limiter
{
	if(x<1)
		return x;
	else
		return 1;		
}

double miu_Ven(double x)//Venkatakrishnan limiter
{
	return (x*x+2*x)/(x*x+x+2);
}
