#include <stdio.h>
#include <stdlib.h>

#define BUF_SIZE 500000

void display( int array[], int maxlen)
{
	int i;
	for( i=0;i<maxlen;i++)
	{
		printf("%-5d",array[i]);
	}
	printf("\n");

	return ;
}


void swap(int *a, int *b)
{
	int temp;
	
	temp = *a;
	*a = *b;
	*b = temp;
	return;
}

void quicksort(int array[], int maxlen, int begin, int end)
{
	int i,j;
	if(begin < end )
	{
		i = begin + 1;
		j = end;
		while( i < j )
		{
			if(array[i]>array[begin])
			{
				swap( &array[i], &array[j]);
				j--;
			}
			else 
			{
				i++;
			}
		}

		if(array[i] >= array[begin] )
		{
			i--;
		}
		swap( &array[begin], &array[i]);

		quicksort(array, maxlen, begin, i);
		quicksort(array, maxlen, j, end);
	}
}



int main()
{
	int n;
	int array[BUF_SIZE];
	int maxlen = BUF_SIZE;
	int i;

	for(i = 0; i < maxlen; i++)
	{
		array[i] = maxlen-i;
	}
	printf("Before sort\n");
	//display(array, maxlen);
	
	quicksort(array, maxlen, 0 , maxlen - 1);

	printf("After sort\n");
	//display(array, maxlen);

	return 0;
	 
}







