void set_receiver(int rec_num, int node_num, double edge_size,  double *rec_x, double *rec_y, double **node_xy, double *rec_node)
/******************************************************************************/
/*
  Purpose:

    set_receiver return receivers node index according to the given (rec_x, rec_y), 
    the total number of receiver is rec_num.

*/
{
	int i;
	double x, y;
	if ( rec_num >=1)
	{
	  for( i = 0; i < rec_num; i++)
	  {	
	        rec_node[i] = -1; // initial to be -1
		x = rec_x[i];
		y = rec_y[j];
		rec_node[i] = node_location( node_xy, x, y, edge_size, node_num );
	  }
	}
	else 
	{
	  printf("Error: receiver number is less than 1.\n");
	}

}
