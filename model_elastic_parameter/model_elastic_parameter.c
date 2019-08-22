
void model_elastic_parameter(int node_num, int element_num, int element_order, int *element_node, double **node_xy, double *rho, double *vp, double *vs, double **c)
/******************************************************************************/
/*
  Purpose:

   elastic_model is used to set the model and gives rho[node_num], vp[node_num], vs[node_num], c[4][node_num]: c11, c13, c33, c44. 
   
   attention: the rho, c11, c12, c33 and c44 can aslo be set in the assblebing process by giving their values according their x_val and y_val location,
   which allows to set diffreent values even in the same element.  Because their valu, es can be defined on the qudrature rule points instead of the element control points.

*/

{

	int p, order, element;
	double x, y;
	double RHO, VP, VS, lamda, miu;

	double radius = 20.0;
	double c_center_x = 1000; // 1000.0 - 60.0;
	double c_center_y = 940.0;

	for (element = 0; element < element_num * element_order; element = element + element_order)
	{

		for (order = 0; order < element_order; order++)
		{
			// p is the current node and set rho, vp, vs according to its coordinates(x,y).

			p = element_node[order + element] - 1; // convert to 0-based index
			x = node_xy[0][p];
			y = node_xy[1][p];
			RHO = 2200.0;
			VP = 2000.0;
			VS = 1154.7;
			miu = VS * VS * RHO;
			lamda = VP * VP * RHO - 2 * miu;
			rho[p] = RHO;
			vp[p] = VP;
			vs[p] = VS;
			c[0][p] = lamda + 2 * miu;
			c[1][p] = lamda;
			c[2][p] = lamda + 2 * miu;
			c[3][p] = miu;
			/*	if ( (x - c_center_x)*(x - c_center_x) + (y - c_center_y)*(y - c_center_y) > radius * radius ) //  set model condition
			{
				RHO = 2200.0;
				VP = 2000.0;
				VS = 1154.7;
				miu = VS * VS * RHO;
				lamda = VP * VP * RHO - 2 * miu;
				rho[p] = RHO;
                vp[p] = VP;
                vs[p] = VS;
				c[0][p] = lamda + 2 * miu;
				c[1][p] = lamda;
				c[2][p] = lamda + 2 * miu;
				c[3][p] = miu;
			}
			else
			{
				RHO = 2200.0;
				VP = 500.0;
				VS = 288.7;
                miu = VS * VS * RHO;
                lamda = VP * VP * RHO - 2 * miu;
                rho[p] = RHO;
                vp[p] = VP;
                vs[p] = VS;
                c[0][p] = lamda + 2 * miu;
                c[1][p] = lamda;
                c[2][p] = lamda + 2 * miu;
                c[3][p] = miu;
			}
			*/
		}
	}
}
