#include <stdio.h>
#include <stdlib.h>
#define n 1000000
int main()


{


	int ai[n]; //={2,1,1,4,5,6};
	int i,j;
	int temp;
	
	
	for(i = 0; i < n; i++)
	{
		ai[i]=n-i;
	//	printf("%d ", ai[i] );
	}
	//printf("\n"); 

	for(i = 0; i < n-1; i++)
	{
		for(j = 0; j < n-1-i; j++)
		{
			if(ai[j] > ai[j+1])
			{
				temp = ai[j];
				ai[j] = ai[j+1];
				ai[j+1] = temp;
			}
		}
	}

	/*for(i = 0; i < n; i++)
	{
		printf("%d ", ai[i] );
	}
	printf("\n");
	*/
}
