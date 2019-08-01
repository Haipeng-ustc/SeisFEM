#include <stdio.h>
#include <stdlib.h>

#define N 100000
int main()
{
        //int ia[N] = {0,0,1,1,2,2,3,3,3};
	//int ja[N] = {0,1,1,2,4,2,3,1,1};
	//double a[N] = {1.0,7.0,2.0,8.0,5.0,3.0,9.0,6.0,4.0};

	int ia[N] = {0};
	int ja[N] = {0};
        int ib[N] = {0};
	int jb[N] = {0};
	double b[N] = {0.0};
	double a[N] = {0.0};
	int flag[N] ={0};
	int i,j;
	
	for(i = 0; i < N; i++)
	{
		ia[i] = 1;
		ja[i] = 1;
		 a[i] = 1.0;
		if(i % 10 == 0 )
		{
		ia[i] = 10;
		ja[i] = 10;
		 a[i] = 2.0;
		}
	}	
	

	printf("Before sum the duplicate element \n");
	/*for ( i = 0; i < N;i++)
	{
		printf("%5d %5d 	%8f %5d \n", ia[i], ja[i], a[i], flag[i]);
 	}
	*/
	for(i = 0;i < N - 1; i++)
	{
		for(j = i + 1; j<N; j++)
		{ 

			if( ia[i] == ia[j] && ja[i] == ja[j] && flag[j] == 0 )
			{
				a[i] = a[i] + a[j];		
				flag[j] = 1;
			}		
		}
  	} 
	
	printf("After sum the duplicate element \n");
	/* for ( i = 0; i < N; i++)
	{
		printf("%5d %5d 	%8f %5d \n", ia[i], ja[i], a[i], flag[i]);
 	}
	*/
	j = 0;

	for(i = 0; i < N; i++)
	{
		if (flag[i] == 0)
		{

			ib[j] = ia[i];
			jb[j] = ja[i];
			 b[j] =  a[i];
			j++;
		}
	} 
	printf("Final Result, total number is %d \n",j);
	for(i = 0; i < j; i++)
	{
		printf("%5d %5d 	%8f\n", ib[i], jb[i], b[i]);
	}

	return 0;
}
