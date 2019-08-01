
void elastic_model(int node_num, int element_num, int element_order, int *element_node, double **node_xy, double *rho, double **c)
/******************************************************************************/
/*
  Purpose:

   elastic_model is used to set the model and gives rho[node_num], c[4][node_num]: c11, c13, c33, c44. 
   
   attention: the rho, c11, c12, c33 and c44 can aslo be set in the assblebing process by giving their values according their x_val and y_val location,
   which allows to set diffreent values even in the same element.  Because their valu, es can be defined on the qudrature rule points instead of the element control points.

*/

{

	int p, order, element;
	double x, y;
	double RHO, vp, vs, lamda, miu;

	for (element = 0; element < element_num * element_order; element = element + element_order)
	{

		for (order = 0; order < element_order; order++)
		{
			// p is the current node and set rho, vp, vs according to its coordinates(x,y).

			p = element_node[order + element] - 1; // convert to 0-based index
			x = node_xy[0][p];
			y = node_xy[1][p];

			if (x >= 0 && y >= 5) //  set model condition
			{
				RHO = 2200;
				vp = 2200;
				vs = 1154.7;
				miu = vs * vs * RHO;
				lamda = vp * vp * RHO - 2 * miu;
				rho[p] = RHO;
				c[0][p] = lamda + 2 * miu;
				c[1][p] = lamda;
				c[2][p] = lamda + 2 * miu;
				c[3][p] = miu;
			}
			else
			{
				RHO = 2200;
				vp = 3000;
				vs = 1732.1;
				miu = vs * vs * RHO;
				lamda = vp * vp * RHO - 2 * miu;
				rho[p] = RHO;
				c[0][p] = lamda + 2 * miu;
				c[1][p] = lamda;
				c[2][p] = lamda + 2 * miu;
				c[3][p] = miu;
			}
		}
	}
}
