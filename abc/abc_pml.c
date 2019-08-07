void abc_pml(int node_num, int element_num, int element_order, int *element_node, double **node_xy,
							int pml_nx, int pml_ny, double h, double xmin, double xmax, double ymin, double ymax, double vp_max,
							int use_pml_xmin, int use_pml_xmax, int use_pml_ymin, int use_pml_ymax,
							double *pml_dx, double *pml_dy, double *pml_dxx, double *pml_dyy)
/******************************************************************************/
/*
  Purpose:

   absorbing_boundary_pml is used to set pml absorbing boundary condition and give pml_dx[node_num], 
   pml_dy[node_num], pml_dxx[node_num], pml_dyy[node_num] for assembing the matrices. 
   
   attention: the  pml_dx[node_num], pml_dy[node_num], pml_dxx[node_num], pml_dyy[node_num] can aslo 
   be set in the assblebing process by giving their values according their x_val and y_val location,
   which allows to set diffreent values even in the same element.  Because their values can be defined
   on the qudrature rule points instead of the element control points.
  
   Output:
     pml_dx, pml_dy, pml_dxx, pml_dyy
   
  
           +--------------------+   ---> x
           | 12 |    2     |23  |
           |--------------------|
           |    |          |    |
           |    |          |    |
           |    |          |    |
           |  1 |          |  3 |
           |    |          |    |
           |    |          |    |
           |    |          |    |
           |    |          |    |
           |--------------------|
           | 14 |    4     | 34 |
           +--------------------+
           |       
           | z
           +
*/
{

	int p, order, element;
	double x, y;
	double d0_x, d0_y;
	double Rcoef;
	double xmin_pml, xmax_pml, ymin_pml, ymax_pml, pml_x_thick, pml_y_thick;

	Rcoef = 0.0001;
	pml_x_thick = h * pml_nx;
	pml_y_thick = h * pml_ny;
	xmin_pml = xmin + pml_x_thick;
	xmax_pml = xmax - pml_x_thick;
	ymin_pml = ymin + pml_y_thick;
	ymax_pml = ymax - pml_y_thick;

	d0_x = 3.0 * vp_max * log10(1.0 / Rcoef) / (2.0 * pml_x_thick);
	d0_y = 3.0 * vp_max * log10(1.0 / Rcoef) / (2.0 * pml_y_thick);

	for (element = 0; element < element_num * element_order; element = element + element_order)
	{

		for (order = 0; order < element_order; order++)
		{
			// p is the current node and set pml_dx, pml_dy, pml_dxx, pml_dyy according to its coordinates(x,y).
			p = element_node[order + element] - 1; // convert to 0-based index
			x = node_xy[0][p];
			y = node_xy[1][p];

			// x direction
			if ((x <= xmin_pml) && (use_pml_xmin == 1))
			{
				pml_dx[p] = d0_x * ((xmin_pml - x) / pml_x_thick) * ((xmin_pml - x) / pml_x_thick);
				pml_dxx[p] = -d0_x * 2 * (xmin_pml - x) / (pml_x_thick * pml_x_thick);
				// Be careful when caculating the derivative, expecially the direction
			}
			else if ((x >= xmax_pml) && (use_pml_xmax == 1))
			{
				pml_dx[p] = d0_x * ((x - xmax_pml) / pml_x_thick) * ((x - xmax_pml) / pml_x_thick);
				pml_dxx[p] = d0_x * 2 * (x - xmax_pml) / (pml_x_thick * pml_x_thick);
			}
			else
			{
				pml_dx[p] = 0.0;
				pml_dxx[p] = 0.0;
			}

			// y direction
			if ((y <= ymin_pml) && (use_pml_ymin == 1))
			{
				pml_dy[p] = d0_y * ((ymin_pml - y) / pml_y_thick) * ((ymin_pml - y) / pml_y_thick);
				pml_dyy[p] = -d0_y * 2 * (ymin_pml - y) / (pml_y_thick * pml_y_thick);
				// Be careful when caculating the derivative, expecially the direction
			}
			else if ((y >= ymax_pml) && (use_pml_ymax == 1))
			{
				pml_dy[p] = d0_y * ((y - ymax_pml) / pml_y_thick) * ((y - ymax_pml) / pml_y_thick);
				pml_dyy[p] = d0_y * 2 * (y - ymax_pml) / (pml_y_thick * pml_y_thick);
			}
			else
			{
				pml_dy[p] = 0.0;
				pml_dyy[p] = 0.0;
			}
		}
	}
}