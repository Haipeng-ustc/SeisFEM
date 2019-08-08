void set_receiver_node(int rec_num, int node_num, double edge_size, double *rec_x, double *rec_y, double **node_xy, int *rec_node)
/******************************************************************************/
/*
  Purpose:

    set_receiver_node return receivers node index according to the given (rec_x, rec_y), 
    the total number of receiver is rec_num.

*/
{
	int i;
	double x, y;
	FILE *fp_receiver_xy;
	fp_receiver_xy = fopen("receiver_xy.dat", "w");
	if (rec_num >= 1)
	{
		for (i = 0; i < rec_num; i++)
		{
			rec_node[i] = -1; // initial to be -1
			x = rec_x[i];
			y = rec_y[i];
			rec_node[i] = node_location(node_num, edge_size, node_xy, x, y);
			fprintf(fp_receiver_xy, "%lf	%lf\n", node_xy[0][rec_node[i]], node_xy[1][rec_node[i]]);
		}
	}
	else
	{
		printf("Error: receiver number is less than 1.\n");
	}
	fclose(fp_receiver_xy);
}
