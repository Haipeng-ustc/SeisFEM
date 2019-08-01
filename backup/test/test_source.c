#include <stdio.h>
#include <math.h>
#include "../source_receiver/seismic_source.c"
int main()
{

	double f0, t0, t, point_source;
	int i;
	f0 = 30; 
	t0 =1.2/f0;
	t  = 0.005;
	for(i = 0; i < 100;i++)

	{
		point_source = seismic_source(f0,t0,1.0,t*i); 
		printf("%f	",point_source);
	}
	return 0;



}
