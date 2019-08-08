void set_source_node(int src_num, int node_num, double edge_size, double *src_x, double *src_y, double **node_xy, int *src_node)
/******************************************************************************/
/*
  Purpose:

    set_source_node return source node index according to the given (src_x, src_y), 
    the total number of source is src_num.

*/
{
	int i;
	double x, y;
	FILE *fp_source_xy;
	fp_source_xy = fopen("source_xy.dat", "w");
	if (src_num >= 1)
	{
		for (i = 0; i < src_num; i++)
		{
			src_node[i] = -1; // initial to be -1
			x = src_x[i];
			y = src_y[i];
			src_node[i] = node_location(node_num, edge_size, node_xy, x, y);
			fprintf(fp_source_xy, "%lf	%lf\n", node_xy[0][src_node[i]], node_xy[1][src_node[i]]);
		}
	}
	else
	{
		printf("Error: source number is less than 1.\n");
	}
	fclose(fp_source_xy);
}
