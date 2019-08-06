
void acoustic_model(int node_num, int element_num, int element_order, int *element_node, double **node_xy, double *rho, double *vp)
/******************************************************************************/
/*
  Purpose:

   acoustic_model is used to set the model and gives rho[node_num], vp[node_num]. 
   
   attention: the rho, vp, vs can aslo be set in the assblebing process by giving
   their values according their x_val and y_val location, which allows to set 
   diffreent values even in the same element.  Because their values can be defined
   on the qudrature rule points instead of the element control points.

*/

{

	int p, order, element;
	double x, y;

	for (element = 0; element < element_num * element_order; element = element + element_order)
	{

		for (order = 0; order < element_order; order++)
		{
			// p is the current node and set rho, vp, vs according to its coordinates(x,y).

			p = element_node[order + element] - 1; // convert to 0-based index
			x = node_xy[0][p];
			y = node_xy[1][p];

			if (x >= 0 && y >= 0) //  set model condition
			{
				rho[p] = 2200;
				vp[p] = 2000;
			}
			else
			{
				rho[p] = 2200;
				vp[p] = 3000;
			}
		}
	}
}
