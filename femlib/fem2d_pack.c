# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "fem2d_pack.h"

/******************************************************************************/

void bandwidth_mesh ( int element_order, int element_num, int element_node[],
  int *ml, int *mu, int *m )

/******************************************************************************/
/*
  Purpose:

    BANDWIDTH_MESH determines the bandwidth of the coefficient matrix.

  Discussion:

    The quantity computed here is the "geometric" bandwidth determined
    by the finite element mesh alone.

    If a single finite element variable is associated with each node
    of the mesh, and if the nodes and variables are numbered in the
    same way, then the geometric bandwidth is the same as the bandwidth
    of a typical finite element matrix.

    The bandwidth M is defined in terms of the lower and upper bandwidths:

      M = ML + 1 + MU

    where 

      ML = maximum distance from any diagonal entry to a nonzero
      entry in the same row, but earlier column,

      MU = maximum distance from any diagonal entry to a nonzero
      entry in the same row, but later column.

    Because the finite element node adjacency relationship is symmetric,
    we are guaranteed that ML = MU.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 January 2006

  Author:

    John Burkardt

  Parameters:

    Input, int ELEMENT_ORDER, the order of the elements.

    Input, int ELEMENT_NUM, the number of elements.

    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM];
    ELEMENT_NODE(I,J) is the global index of local node I in element J.

    Output, int *ML, *MU, the lower and upper bandwidths of the matrix.

    Output, int *M, the bandwidth of the matrix.
*/
{
  int element;
  int global_i;
  int global_j;
  int local_i;
  int local_j;

  *ml = 0;
  *mu = 0;

  for ( element = 0; element < element_num; element++ )
  {
    for ( local_i = 0; local_i < element_order; local_i++ )
    {
      global_i = element_node[local_i+element*element_order];

      for ( local_j = 0; local_j < element_order; local_j++ )
      {
        global_j = element_node[local_j+element*element_order];

        *mu = i4_max ( *mu, global_j - global_i );
        *ml = i4_max ( *ml, global_i - global_j );
      }
    }
  }

  *m = *ml + 1 + *mu;

  return;
}
/******************************************************************************/

void bandwidth_var ( int element_order, int element_num, int element_node[],
  int node_num, int var_node[], int var_num, int var[], int *ml, int *mu, 
  int *m )

/******************************************************************************/
/*
  Purpose:

    BANDWIDTH_VAR determines the bandwidth for finite element variables.

  Discussion:

    We assume that, attached to each node in the finite element mesh
    there are a (possibly zero) number of finite element variables.
    We wish to determine the bandwidth necessary to store the stiffness
    matrix associated with these variables.

    An entry K(I,J) of the stiffness matrix must be zero unless the
    variables I and J correspond to nodes N(I) and N(J) which are
    common to some element.

    In order to determine the bandwidth of the stiffness matrix, we
    essentially seek a nonzero entry K(I,J) for which abs ( I - J )
    is maximized.

    The bandwidth M is defined in terms of the lower and upper bandwidths:

      M = ML + 1 + MU

    where

      ML = maximum distance from any diagonal entry to a nonzero
      entry in the same row, but earlier column,

      MU = maximum distance from any diagonal entry to a nonzero
      entry in the same row, but later column.

    We assume the finite element variable adjacency relationship is 
    symmetric, so we are guaranteed that ML = MU.

    Note that the user is free to number the variables in any way
    whatsoever, and to associate variables to nodes in any way,
    so that some nodes have no variables, some have one, and some
    have several.  

    The storage of the indices of the variables is fairly simple.
    In VAR, simply list all the variables associated with node 1, 
    then all those associated with node 2, and so on.  Then set up
    the pointer array VAR_NODE so that we can jump to the section of
    VAR where the list begins for any particular node.

    The routine does not check that each variable is only associated
    with a single node.  This would normally be the case in a finite
    element setting.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 September 2006

  Author:

    John Burkardt

  Parameters:

    Input, int ELEMENT_ORDER, the order of the elements.

    Input, int ELEMENT_NUM, the number of elements.

    Input,  ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM];
    ELEMENT_NODE(I,J) is the global index of local node I in element J.

    Output, int *ML, *MU, the lower and upper bandwidths of the matrix.

    Output, int *M, the bandwidth of the matrix.
*/
{
  int element;
  int node_global_i;
  int node_global_j;
  int node_local_i;
  int node_local_j;
  int var_global_i;
  int var_global_j;
  int var_local_i;
  int var_local_j;

  *ml = 0;
  *mu = 0;

  for ( element = 0; element < element_num; element++ )
  {
    for ( node_local_i = 0; node_local_i < element_order; node_local_i++ )
    {
      node_global_i = element_node[node_local_i+element*element_order];

      for ( var_local_i = var_node[node_global_i-1]; 
            var_local_i <= var_node[node_global_i]-1; var_local_i++ )
      {
        var_global_i = var[var_local_i-1];

        for ( node_local_j = 0; node_local_j < element_order; node_local_j++ )
        {
          node_global_j = element_node[node_local_j+element*element_order];

          for ( var_local_j = var_node[node_global_j-1]; 
                var_local_j <= var_node[node_global_j]-1; var_local_j++ )
          {
            var_global_j = var[var_local_j-1];

            *mu = i4_max ( *mu, var_global_j - var_global_i );
            *ml = i4_max ( *ml, var_global_i - var_global_j );
          }
        }
      }
    }
  }

  *m = *ml + 1 + *mu;

  return;
}
/******************************************************************************/

void basis_11_t3 ( double t[2*3], int i, double p[2], double *qi, 
  double *dqidx, double *dqidy )

/******************************************************************************/
/*
  Purpose:

    BASIS_11_T3: one basis at one point for a T3 element.

  Discussion:

    The routine is given the coordinates of the nodes of a triangle.

           3
          / \
         /   \
        /     \
       1-------2

    It evaluates the linear basis function Q(I)(X,Y) associated with
    node I, which has the property that it is a linear function
    which is 1 at node I and zero at the other two nodes.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 January 2006

  Author:

    John Burkardt

  Parameters:

    Input, double T[2*3], the coordinates of the nodes.

    Input, int I, the index of the desired basis function.
    I should be between 1 and 3.

    Input, double P[2], the coordinates of the point where 
    the basis function is to be evaluated.

    Output, double *QI, *DQIDX, *DQIDY, the value of the I-th basis function
    and its X and Y derivatives.
*/
{
  double area;
  int ip1;
  int ip2;

  area = t[0+0*2] * ( t[1+1*2] - t[1+2*2] ) 
       + t[0+1*2] * ( t[1+2*2] - t[1+0*2] ) 
       + t[0+2*2] * ( t[1+0*2] - t[1+1*2] );

  if ( area == 0.0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "BASIS_11_T3 - Fatal error!\n" );
    fprintf ( stderr, "  Element has zero area.\n" );
    fprintf ( stderr, "  Area = %g\n", area );
    exit ( 1 );
  }

  if ( i < 1 || 3 < i )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "BASIS_11_T3 - Fatal error!\n" );
    fprintf ( stderr, "  Basis index I is not between 1 and 3.\n" );
    fprintf ( stderr, "  I = %d\n", i );
    exit ( 1 );
  }

  ip1 = i4_wrap ( i + 1, 1, 3 );
  ip2 = i4_wrap ( i + 2, 1, 3 );

  *qi = ( ( t[0+(ip2-1)*2] - t[0+(ip1-1)*2] ) 
        * ( p[1]           - t[1+(ip1-1)*2] ) 
        - ( t[1+(ip2-1)*2] - t[1+(ip1-1)*2] ) 
        * ( p[0]           - t[0+(ip1-1)*2] ) ) / area;

  *dqidx = - ( t[1+(ip2-1)*2] - t[1+(ip1-1)*2] ) / area;
  *dqidy =   ( t[0+(ip2-1)*2] - t[0+(ip1-1)*2] ) / area;

  return;
}
/******************************************************************************/

void basis_11_t3_test ( )

/******************************************************************************/
/*
  Purpose:

    BASIS_11_T3_TEST verifies BASIS_11_T3.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 January 2006

  Author:

    John Burkardt

  Parameters:

    None
*/
{
# define NODE_NUM 3

  double dqjdx;
  double dqjdy;
  int i;
  int j;
  double qj;
  double p[2];
  double sum_x;
  double sum_y;
  double t[2*NODE_NUM] = { 
    2.0, 0.0,
    4.0, 3.0,
    0.0, 4.0 };

  printf ( "\n" );
  printf ( "BASIS_11_T3_TEST:\n" );
  printf ( "  Verify basis functions for element T3.\n" );
  printf ( "\n" );
  printf ( "  Number of nodes = %d\n", NODE_NUM );

  printf ( "\n" );
  printf ( "  Physical Nodes:\n" );
  printf ( "\n" );
  for ( j = 0; j < NODE_NUM; j++ )
  {
    printf ( "  %3d  %10g  %10g\n", j, t[0+j*2], t[1+j*2] );
  }
 
  printf ( "\n" );
  printf ( "  The basis function values at basis nodes\n" );
  printf ( "  should form the identity matrix.\n" );
  printf ( "\n" );

  for ( i = 0; i < NODE_NUM; i++ )
  {
    p[0] = t[0+i*2];
    p[1] = t[1+i*2];

    for ( j = 0; j < NODE_NUM; j++ )
    {
      basis_11_t3 ( t, j+1, p, &qj, &dqjdx, &dqjdy );
      printf ( "  %10g", qj );
    }
    printf ( "\n" );
  }

  printf ( "\n" );
  printf ( "  The X and Y derivatives should sum to 0.\n" );
  printf ( "\n" );
  printf ( "  dPhidX sum, dPhidY sum:\n" );
  printf ( "\n" );

  for ( i = 0; i < NODE_NUM; i++ )
  {
    p[0] = t[0+i*2];
    p[1] = t[1+i*2];

    sum_x = 0.0;
    sum_y = 0.0;
    for ( j = 0; j < NODE_NUM; j++ )
    {
      basis_11_t3 ( t, j+1, p, &qj, &dqjdx, &dqjdy );
      sum_x = sum_x + dqjdx;
      sum_y = sum_y + dqjdy;
    }
    printf ( "  %10g  %10g\n", sum_x, sum_y );
  }

  return;
# undef NODE_NUM
}
/******************************************************************************/

void basis_11_t4 ( double t[2*4], int i, double p[], double *phi, 
  double *dphidx, double *dphidy )

/******************************************************************************/
/*
  Purpose:

    BASIS_MN_T4: one basis at one point for a T4 element.

  Discussion:

    The T4 element is the cubic bubble triangle.

    The routine is given the coordinates of the vertices of a triangle.
    It works directly with these coordinates, and does not refer to a 
    reference element.

    The sides of the triangle DO NOT have to lie along a coordinate
    axis.

    The routine evaluates the basis functions associated with each vertex,
    and their derivatives with respect to X and Y.

  Physical Element T4: 
       
            3
           / \
          /   \
         /  4  \
        /       \
       1---------2

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 August 2009

  Author:

    John Burkardt

  Parameters:

    Input, double T[2*4], the coordinates of the vertices
    of the triangle, and the coordinates of the centroid.  
    It is common to list the first three points in counter clockwise
    order.

    Input, int I, the index of the basis function.

    Input, double P[2], the point where the basis function
    is to be evaluated.

    Output, double *PHI, the value of the basis function
    at the evaluation point.

    Output, double *DPHIDX, *DPHIDY, the value of the 
    derivatives at the evaluation point.

  Local parameters:

    Local, double AREA, is (twice) the area of the triangle.
*/
{
  double area;
  double dpsidx[4];
  double dpsidy[4];
  int j;
  double psi[4];

  area = t[0+0*2] * ( t[1+1*2] - t[1+2*2] ) 
       + t[0+1*2] * ( t[1+2*2] - t[1+0*2] ) 
       + t[0+2*2] * ( t[1+0*2] - t[1+1*2] );

  psi[0] =     (   ( t[0+2*2] - t[0+1*2] ) * ( p[1] - t[1+1*2] )     
                 - ( t[1+2*2] - t[1+1*2] ) * ( p[0] - t[0+1*2] ) );
  dpsidx[0] =    - ( t[1+2*2] - t[1+1*2] );
  dpsidy[0] =      ( t[0+2*2] - t[0+1*2] );

  psi[1] =     (   ( t[0+0*2] - t[0+2*2] ) * ( p[1] - t[1+2*2] )     
                 - ( t[1+0*2] - t[1+2*2] ) * ( p[0] - t[0+2*2] ) );
  dpsidx[1] =    - ( t[1+0*2] - t[1+2*2] );
  dpsidy[1] =      ( t[0+0*2] - t[0+2*2] );

  psi[2] =     (   ( t[0+1*2] - t[0+0*2] ) * ( p[1] - t[1+0*2] )     
                 - ( t[1+1*2] - t[1+0*2] ) * ( p[0] - t[0+0*2] ) );
  dpsidx[2] =    - ( t[1+1*2] - t[1+0*2] );
  dpsidy[2] =      ( t[0+1*2] - t[0+0*2] );
/*
  Normalize the first three functions.
*/
  for ( j = 0; j < 3; j++ )
  {
    psi[j]    = psi[j]    / area;
    dpsidx[j] = dpsidx[j] / area;
    dpsidy[j] = dpsidy[j] / area;
  }
/*
  Compute the cubic bubble function.
*/
  psi[3] = 27.0 * psi[0] * psi[1] * psi[2];

  dpsidx[3] = 27.0 * (
                  dpsidx[0] *    psi[1] *    psi[2]
                  +  psi[0] * dpsidx[1] *    psi[2]
                  +  psi[0] *    psi[1] * dpsidx[2] );

  dpsidy[3] = 27.0 * (
                  dpsidy[0] *    psi[1] *    psi[2]
                  +  psi[0] * dpsidy[1] *    psi[2]
                  +  psi[0] *    psi[1] * dpsidy[2] );
/*
  Subtract 1/3 of the cubic bubble function from each of the three linears.
*/
  for ( j = 0; j < 3; j++ )
  {
       psi[j] =    psi[j] -    psi[3] / 3.0;
    dpsidx[j] = dpsidx[j] - dpsidx[3] / 3.0;
    dpsidy[j] = dpsidy[j] - dpsidy[3] / 3.0;
  }

  *phi = psi[i-1];
  *dphidx = dpsidx[i-1];
  *dphidy = dpsidy[i-1];

  return;
}
/******************************************************************************/

void basis_11_t4_test ( )

/******************************************************************************/
/*
  Purpose:

    BASIS_11_T4_TEST verifies BASIS_11_T4.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 March 2009

  Author:

    John Burkardt

  Parameters:

    None
*/
{
# define NODE_NUM 4

  double dqjdx;
  double dqjdy;
  int i;
  int j;
  double qj;
  double p[2];
  double sum_x;
  double sum_y;
  double t[2*NODE_NUM] = { 
    2.0, 0.0,
    4.0, 3.0,
    0.0, 4.0,
    0.0, 0.0 };
/*
  The node associated with the fourth basis function is the centroid.
*/
  t[0+3*2] = ( t[0+0*2] + t[0+1*2] + t[0+2*2] ) / 3.0;
  t[1+3*2] = ( t[1+0*2] + t[1+1*2] + t[1+2*2] ) / 3.0;

  printf ( "\n" );
  printf ( "BASIS_11_T4_TEST:\n" );
  printf ( "  Verify basis functions for element T4.\n" );
  printf ( "\n" );
  printf ( "  Number of nodes = %d\n", NODE_NUM );

  printf ( "\n" );
  printf ( "  Physical Nodes:\n" );
  printf ( "\n" );
  for ( j = 0; j < NODE_NUM; j++ )
  {
    printf ( "  %3d  %10g  %10g\n", j, t[0+j*2], t[1+j*2] );
  }
 
  printf ( "\n" );
  printf ( "  The basis function values at basis nodes\n" );
  printf ( "  should form the identity matrix.\n" );
  printf ( "\n" );

  for ( i = 0; i < NODE_NUM; i++ )
  {
    p[0] = t[0+i*2];
    p[1] = t[1+i*2];

    for ( j = 0; j < NODE_NUM; j++ )
    {
      basis_11_t4 ( t, j+1, p, &qj, &dqjdx, &dqjdy );
      printf ( "  %10g", qj );
    }
    printf ( "\n" );
  }

  printf ( "\n" );
  printf ( "  The X and Y derivatives should sum to 0.\n" );
  printf ( "\n" );
  printf ( "  dPhidX sum, dPhidY sum:\n" );
  printf ( "\n" );

  for ( i = 0; i < NODE_NUM; i++ )
  {
    p[0] = t[0+i*2];
    p[1] = t[1+i*2];

    sum_x = 0.0;
    sum_y = 0.0;
    for ( j = 0; j < NODE_NUM; j++ )
    {
      basis_11_t4 ( t, j+1, p, &qj, &dqjdx, &dqjdy );
      sum_x = sum_x + dqjdx;
      sum_y = sum_y + dqjdy;
    }
    printf ( "  %10g  %10g\n", sum_x, sum_y );
  }

  return;
# undef NODE_NUM
}
/******************************************************************************/

void basis_11_t6 ( double t[2*6], int i, double p[], double *bi, 
  double *dbidx, double *dbidy )

/******************************************************************************/
/*
  Purpose:

    BASIS_11_T6: one basis at one point for the T6 element.

  Discussion:

    The routine is given the coordinates of the nodes of a triangle. 
       
           3
          / \
         6   5
        /     \
       1---4---2

    It evaluates the quadratic basis function B(I)(X,Y) associated with
    node I, which has the property that it is a quadratic function
    which is 1 at node I and zero at the other five nodes.

    This routine assumes that the sides of the triangle are straight,
    so that the midside nodes fall on the line between two vertices.

    This routine relies on the fact that each basis function can be
    written as the product of two linear factors, which are easily
    computed and normalized.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    02 February 2006

  Author:

    John Burkardt

  Parameters:

    Input, double T[2*6], the coordinates of the nodes.

    Input, int I, the index of the desired basis function.
    I should be between 1 and 6.

    Input, double P[2], the coordinates of a point at which the basis
    function is to be evaluated.

    Output, double *BI, *DBIDX, *DBIDY, the values of the basis function
    and its X and Y derivatives.
*/
{
  double gf;
  double gn;
  double hf;
  double hn;
  int j1;
  int j2;
  int k1;
  int k2;

  if ( i < 1 || 6 < i )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "BASIS_11_T6 - Fatal error!\n" );
    fprintf ( stderr, "  Basis index I is not between 1 and 6.\n" );
    fprintf ( stderr, "  I = %d\n", i );
    exit ( 1 );
  }
/*
  Determine the pairs of nodes.
*/
  if ( i <= 3 )
  {
    j1 = i4_wrap ( i + 1, 1, 3 );
    j2 = i4_wrap ( i + 2, 1, 3 );
    k1 = i + 3;
    k2 = i4_wrap ( i + 5, 4, 6 );
  }
  else
  {
    j1 = i - 3;
    j2 = i4_wrap ( i - 3 + 2, 1, 3 );
    k1 = i4_wrap ( i - 3 + 1, 1, 3 );
    k2 = i4_wrap ( i - 3 + 2, 1, 3 );
  }
/*
  For C indexing, it is helpful to knock the indices down by one.
*/
  i  = i  - 1;
  j1 = j1 - 1;
  j2 = j2 - 1;
  k1 = k1 - 1;
  k2 = k2 - 1;
/*
  Evaluate the two linear factors GF and HF, 
  and their normalizers GN and HN.
*/
  gf = ( p[0]      - t[0+j1*2] ) * ( t[1+j2*2] - t[1+j1*2] ) 
     - ( t[0+j2*2] - t[0+j1*2] ) * ( p[1]      - t[1+j1*2] );

  gn = ( t[0+i*2]  - t[0+j1*2] ) * ( t[1+j2*2] - t[1+j1*2] ) 
     - ( t[0+j2*2] - t[0+j1*2] ) * ( t[1+i*2]  - t[1+j1*2] );

  hf = ( p[0]      - t[0+k1*2] ) * ( t[1+k2*2] - t[1+k1*2] ) 
     - ( t[0+k2*2] - t[0+k1*2] ) * ( p[1]      - t[1+k1*2] );

  hn = ( t[0+i*2]  - t[0+k1*2] ) * ( t[1+k2*2] - t[1+k1*2] ) 
     - ( t[0+k2*2] - t[0+k1*2] ) * ( t[1+i*2]  - t[1+k1*2] );
/*
  Construct the basis function and its derivatives.
*/
  *bi =        ( gf                      / gn ) 
             * ( hf                      / hn );

  *dbidx =   ( ( t[1+j2*2] - t[1+j1*2] ) / gn ) 
           * (   hf                      / hn )
           + (   gf                      / gn ) 
           * ( ( t[1+k2*2] - t[1+k1*2] ) / hn );

  *dbidy = - ( ( t[0+j2*2] - t[0+j1*2] ) / gn ) 
           * (   hf                      / hn )
           - (   gf                      / gn ) 
           * ( ( t[0+k2*2] - t[0+k1*2] ) / hn );

  return;
}
/******************************************************************************/

void basis_11_t6_test ( )

/******************************************************************************/
/*
  Purpose:

    BASIS_11_T6_TEST verifies BASIS_11_T6.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 February 2006

  Author:

    John Burkardt

  Parameters:

    None
*/
{
# define NODE_NUM 6

  double dphidx[NODE_NUM*NODE_NUM];
  double dphidy[NODE_NUM*NODE_NUM];
  int i;
  int j;
  double p[2];
  double phi[NODE_NUM*NODE_NUM];
  double sum_x;
  double sum_y;
  double t[2*NODE_NUM] = {
    2.0, 0.0,
    4.0, 3.0,
    0.0, 4.0,
    3.0, 1.5,
    2.0, 3.5,
    1.0, 2.0 };
  double v1;
  double v2;
  double v3;

  printf ( "\n" );
  printf ( "BASIS_11_T6_TEST:\n" );
  printf ( "  Verify basis functions for element T6.\n" );
  printf ( "\n" );
  printf ( "  Number of nodes = %d\n", NODE_NUM );

  printf ( "\n" );
  printf ( "  Physical Nodes:\n" );
  printf ( "\n" );
  for ( j = 0; j < NODE_NUM; j++ )
  {
    printf ( "  %3d  %10g  %10g\n", j, t[0+j*2], t[1+j*2] );
  }
 
  printf ( "\n" );
  printf ( "  The basis function values at basis nodes\n" );
  printf ( "  should form the identity matrix.\n" );
  printf ( "\n" );
  for ( i = 1; i <= NODE_NUM; i++ )
  {
    for ( j = 0; j < NODE_NUM; j++ )
    {
      p[0] = t[0+j*2];
      p[1] = t[1+j*2];

      basis_11_t6 ( t, i, p, &v1, &v2, &v3 );

      phi[i-1+j*NODE_NUM] = v1;
      dphidx[i-1+j*NODE_NUM] = v2;
      dphidy[i-1+j*NODE_NUM] = v3;
    }
  }
  for ( i = 0; i < NODE_NUM; i++ )
  {
    for ( j = 0; j < NODE_NUM; j++ )
    {
      printf ( "  %7g", phi[i+j*NODE_NUM] );
    }
    printf ( "\n" );
  }

  printf ( "\n" );
  printf ( "  The X and Y derivatives should sum to 0.\n" );
  printf ( "\n" );
  printf ( "  dPhidX sum, dPhidY sum:\n" );
  printf ( "\n" );
  for ( j = 0; j < NODE_NUM; j++ )
  {
    sum_x = 0.0;
    for ( i = 0; i < NODE_NUM; i++ )
    {
      sum_x = sum_x + dphidx[i+j*NODE_NUM];
    }
    sum_y = 0.0;
    for ( i = 0; i < NODE_NUM; i++ )
    {
      sum_y = sum_y + dphidy[i+j*NODE_NUM];
    }
    printf ( "  %10g  %10g\n", sum_x, sum_y );
  }

  return;
# undef NODE_NUM
}
/******************************************************************************/

void basis_mn_q4 ( double q[2*4], int n, double p[], double phi[], 
  double dphidx[], double dphidy[] )

/******************************************************************************/
/*
  Purpose:

    BASIS_MN_Q4: all bases at N points for a Q4 element.

  Discussion:

    The routine is given the coordinates of the vertices of a quadrilateral.
    It works directly with these coordinates, and does not refer to a 
    reference element.

    The sides of the element are presumed to lie along coordinate axes.

    The routine evaluates the basis functions associated with each corner,
    and their derivatives with respect to X and Y.

  Physical Element Q4:

    |
    |  4-----3
    |  |     |
    |  |     |
    Y  |     |
    |  |     |
    |  |     |
    |  1-----2
    |
    +-----X------>

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 February 2006

  Author:

    John Burkardt

  Parameters:

    Input, double Q[2*4], the coordinates of the vertices.
    It is common to list these points in counter clockwise order.

    Input, int N, the number of evaluation points.

    Input, double P[2*N], the evaluation points.

    Output, double PHI[4*N], the bases at the evaluation points.

    Output, double DPHIDX[4*N], DPHIDY[4*N], the derivatives of the
    bases at the evaluation points.
*/
{
  double area;
  int i;
  int j;

  area =        ( q[0+2*2]         - q[0+0*2] ) 
              * ( q[1+2*2]         - q[1+0*2] );

  for ( j = 0; j < n; j++ )
  {
    phi[0+j*4] =      ( q[0+2*2] - p[0+j*2]             ) 
                    * ( q[1+2*2] - p[1+j*2]             );
    phi[1+j*4] =      (            p[0+j*2]  - q[0+0*2] ) 
                    * ( q[1+2*2] - p[1+j*2]             );
    phi[2+j*4] =      (            p[0+j*2]  - q[0+0*2] )
                    * (            p[1+j*2]  - q[1+0*2] );
    phi[3+j*4] =      ( q[0+2*2] - p[0+j*2]             )
                    * (            p[1+j*2]  - q[1+0*2] );

    dphidx[0+j*4] = - ( q[1+2*2] - p[1+j*2]             );
    dphidx[1+j*4] =   ( q[1+2*2] - p[1+j*2]             );
    dphidx[2+j*4] =   (            p[1+j*2]  - q[1+0*2] );
    dphidx[3+j*4] = - (            p[1+j*2]  - q[1+0*2] );

    dphidy[0+j*4] = - ( q[0+2*2] - p[0+j*2]             );
    dphidy[1+j*4] = - (            p[0+j*2]  - q[0+0*2] );
    dphidy[2+j*4] =   (            p[0+j*2]  - q[0+0*2] );
    dphidy[3+j*4] =   ( q[0+2*2] - p[0+j*2]             );
  }
/*
  Normalize.
*/
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < 4; i++ )
    {
      phi[i+j*4]    = phi[i+j*4]    / area;
      dphidx[i+j*4] = dphidx[i+j*4] / area;
      dphidy[i+j*4] = dphidy[i+j*4] / area;
    }
  }
  return;
}
/******************************************************************************/

void basis_mn_q4_test ( )

/******************************************************************************/
/*
  Purpose:

    BASIS_MN_Q4_TEST verifies BASIS_MN_Q4.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 February 2006

  Author:

    John Burkardt

  Parameters:

    None
*/
{
# define NODE_NUM 4

  double dphidx[NODE_NUM*NODE_NUM];
  double dphidy[NODE_NUM*NODE_NUM];
  int i;
  int j;
  double phi[NODE_NUM*NODE_NUM];
  double q[2*NODE_NUM] = {
    3.0, 1.0,
    5.0, 1.0,
    5.0, 4.0, 
    3.0, 4.0 };
  double sum_x;
  double sum_y;

  printf ( "\n" );
  printf ( "BASIS_MN_Q4_TEST:\n" );
  printf ( "  Verify basis functions for element Q4.\n" );
  printf ( "\n" );
  printf ( "  Number of nodes = %d\n", NODE_NUM );

  printf ( "\n" );
  printf ( "  Physical Nodes:\n" );
  printf ( "\n" );
  for ( i = 0; i < NODE_NUM; i++ )
  {
    printf ( "  %3d  %10g  %10g\n", i, q[0+i*2], q[1+i*2] );
  }
 
  basis_mn_q4 ( q, NODE_NUM, q, phi, dphidx, dphidy );

  printf ( "\n" );
  printf ( "  The basis function values at basis nodes\n" );
  printf ( "  should form the identity matrix.\n" );
  printf ( "\n" );

  for ( i = 0; i < NODE_NUM; i++ )
  {
    for ( j = 0; j < NODE_NUM; j++ )
    {
    printf ( "  %10g", phi[i+j*NODE_NUM] );
    }
    printf ( "\n" );
  }

  printf ( "\n" );
  printf ( "  The X and Y derivatives should sum to 0.\n" );
  printf ( "\n" );
  printf ( "  dPhidX sum, dPhidY sum:\n" );
  printf ( "\n" );

  for ( j = 0; j < NODE_NUM; j++ )
  {
    sum_x = 0.0;
    for ( i = 0; i < NODE_NUM; i++ )
    {
      sum_x = sum_x + dphidx[i+j*NODE_NUM];
    }
    sum_y = 0.0;
    for ( i = 0; i < NODE_NUM; i++ )
    {
      sum_y = sum_y + dphidy[i+j*NODE_NUM];
    }
    printf ( "  %10g  %10g\n", sum_x, sum_y );
  }

  return;
# undef NODE_NUM
}
/******************************************************************************/

void basis_mn_t3 ( double t[2*3], int n, double p[], double phi[], 
  double dphidx[], double dphidy[] )

/******************************************************************************/
/*
  Purpose:

    BASIS_MN_T3: all bases at N points for a T3 element.

  Discussion:

    The routine is given the coordinates of the vertices of a triangle.
    It works directly with these coordinates, and does not refer to a 
    reference element.

    The sides of the triangle DO NOT have to lie along a coordinate
    axis.

    The routine evaluates the basis functions associated with each vertex,
    and their derivatives with respect to X and Y.

  Physical Element T3: 
       
            3
           / \
          /   \
         /     \
        /       \
       1---------2

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 February 2006

  Author:

    John Burkardt

  Parameters:

    Input, double T[2*3], the coordinates of the vertices
    of the triangle.  It is common to list these points in counter clockwise
    order.

    Input, int N, the number of evaluation points.

    Input, double P[2*N], the points where the basis functions 
    are to be evaluated.

    Output, double PHI[3*N], the value of the basis functions 
    at the evaluation points.

    Output, double DPHIDX[3*N], DPHIDY[3*N], the value of the 
    derivatives at the evaluation points.

  Local parameters:

    Local, double AREA, is (twice) the area of the triangle.
*/
{
  double area;
  int i;
  int j;

  area = t[0+0*2] * ( t[1+1*2] - t[1+2*2] ) 
       + t[0+1*2] * ( t[1+2*2] - t[1+0*2] ) 
       + t[0+2*2] * ( t[1+0*2] - t[1+1*2] );

  if ( area == 0.0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "BASIS_MN_T3 - Fatal error!\n" );
    fprintf ( stderr, "  Element has zero area.\n" );
    exit ( 1 );
  }

  for ( j = 0; j < n; j++ )
  {
    phi[0+j*3] =     (   ( t[0+2*2] - t[0+1*2] ) * ( p[1+j*2] - t[1+1*2] )     
                       - ( t[1+2*2] - t[1+1*2] ) * ( p[0+j*2] - t[0+1*2] ) );
    dphidx[0+j*3] =    - ( t[1+2*2] - t[1+1*2] );
    dphidy[0+j*3] =      ( t[0+2*2] - t[0+1*2] );

    phi[1+j*3] =     (   ( t[0+0*2] - t[0+2*2] ) * ( p[1+j*2] - t[1+2*2] )     
                       - ( t[1+0*2] - t[1+2*2] ) * ( p[0+j*2] - t[0+2*2] ) );
    dphidx[1+j*3] =    - ( t[1+0*2] - t[1+2*2] );
    dphidy[1+j*3] =      ( t[0+0*2] - t[0+2*2] );

    phi[2+j*3] =     (   ( t[0+1*2] - t[0+0*2] ) * ( p[1+j*2] - t[1+0*2] )     
                       - ( t[1+1*2] - t[1+0*2] ) * ( p[0+j*2] - t[0+0*2] ) );
    dphidx[2+j*3] =    - ( t[1+1*2] - t[1+0*2] );
    dphidy[2+j*3] =      ( t[0+1*2] - t[0+0*2] );
  }
/*
  Normalize.
*/
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      phi[i+j*3]    = phi[i+j*3]    / area;
      dphidx[i+j*3] = dphidx[i+j*3] / area;
      dphidy[i+j*3] = dphidy[i+j*3] / area;
    }
  }

  return;
}
/******************************************************************************/

void basis_mn_t3_test ( )

/******************************************************************************/
/*
  Purpose:

    BASIS_MN_T3_TEST verifies BASIS_MN_T3.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 February 2006

  Author:

    John Burkardt

  Parameters:

    None.
*/
{
# define NODE_NUM 3

  double dphidx[NODE_NUM*NODE_NUM];
  double dphidy[NODE_NUM*NODE_NUM];
  int i;
  int j;
  double phi[NODE_NUM*NODE_NUM];
  double sum_x;
  double sum_y;
  double t[2*NODE_NUM] = { 
    2.0, 0.0,
    4.0, 3.0,
    0.0, 4.0 };

  printf ( "\n" );
  printf ( "BASIS_MN_T3_TEST:\n" );
  printf ( "  Verify basis functions for element T3.\n" );
  printf ( "\n" );
  printf ( "  Number of nodes = %d\n", NODE_NUM );

  printf ( "\n" );
  printf ( "  Physical Nodes:\n" );
  printf ( "\n" );
  for ( j = 0; j < NODE_NUM; j++ )
  {
    printf ( "  %3d  %10g  %10g\n", j, t[0+j*2], t[1+j*2] );
  }
 
  printf ( "\n" );
  printf ( "  The basis function values at basis nodes\n" );
  printf ( "  should form the identity matrix.\n" );
  printf ( "\n" );

  basis_mn_t3 ( t, NODE_NUM, t, phi, dphidx, dphidy );

  for ( j = 0; j < NODE_NUM; j++ )
  {
    for ( i = 0; i < NODE_NUM; i++ )
    {
      printf ( "  %10g", phi[i+j*NODE_NUM] );
    }
    printf ( "\n" );
  }

  printf ( "\n" );
  printf ( "  The X and Y derivatives should sum to 0.\n" );
  printf ( "\n" );
  printf ( "  dPhidX sum, dPhidY sum:\n" );
  printf ( "\n" );

  for ( j = 0; j < NODE_NUM; j++ )
  {
    sum_x = 0.0;
    sum_y = 0.0;
    for ( i = 0; i < NODE_NUM; i++ )
    {
      sum_x = sum_x + dphidx[i+j*NODE_NUM];
      sum_y = sum_y + dphidy[i+j*NODE_NUM];
    }
    printf ( "  %10g  %10g\n", sum_x, sum_y );
  }

  return;
# undef NODE_NUM
}
/******************************************************************************/

void basis_mn_t4 ( double t[2*4], int n, double p[], double phi[], 
  double dphidx[], double dphidy[] )

/******************************************************************************/
/*
  Purpose:

    BASIS_MN_T4: all bases at N points for a T4 element.

  Discussion:

    The T4 element is the cubic bubble triangle.

    The routine is given the coordinates of the vertices of a triangle.
    It works directly with these coordinates, and does not refer to a 
    reference element.

    The sides of the triangle DO NOT have to lie along a coordinate
    axis.

    The routine evaluates the basis functions associated with each vertex,
    and their derivatives with respect to X and Y.

  Physical Element T4: 
       
            3
           / \
          /   \
         /  4  \
        /       \
       1---------2

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 February 2006

  Author:

    John Burkardt

  Parameters:

    Input, double T[2*4], the coordinates of the vertices
    of the triangle, and the coordinates of the centroid.  
    It is common to list the first three points in counter clockwise
    order.

    Input, int N, the number of evaluation points.

    Input, double P[2*N], the points where the basis functions 
    are to be evaluated.

    Output, double PHI[4*N], the value of the basis functions 
    at the evaluation points.

    Output, double DPHIDX[4*N], DPHIDY[4*N], the value of the 
    derivatives at the evaluation points.

  Local parameters:

    Local, double AREA, is (twice) the area of the triangle.
*/
{
  double area;
  int i;
  int j;

  area = t[0+0*2] * ( t[1+1*2] - t[1+2*2] ) 
       + t[0+1*2] * ( t[1+2*2] - t[1+0*2] ) 
       + t[0+2*2] * ( t[1+0*2] - t[1+1*2] );

  for ( j = 0; j < n; j++ )
  {
    phi[0+j*4] =     (   ( t[0+2*2] - t[0+1*2] ) * ( p[1+j*2] - t[1+1*2] )     
                       - ( t[1+2*2] - t[1+1*2] ) * ( p[0+j*2] - t[0+1*2] ) );
    dphidx[0+j*4] =    - ( t[1+2*2] - t[1+1*2] );
    dphidy[0+j*4] =      ( t[0+2*2] - t[0+1*2] );

    phi[1+j*4] =     (   ( t[0+0*2] - t[0+2*2] ) * ( p[1+j*2] - t[1+2*2] )     
                       - ( t[1+0*2] - t[1+2*2] ) * ( p[0+j*2] - t[0+2*2] ) );
    dphidx[1+j*4] =    - ( t[1+0*2] - t[1+2*2] );
    dphidy[1+j*4] =      ( t[0+0*2] - t[0+2*2] );

    phi[2+j*4] =     (   ( t[0+1*2] - t[0+0*2] ) * ( p[1+j*2] - t[1+0*2] )     
                       - ( t[1+1*2] - t[1+0*2] ) * ( p[0+j*2] - t[0+0*2] ) );
    dphidx[2+j*4] =    - ( t[1+1*2] - t[1+0*2] );
    dphidy[2+j*4] =      ( t[0+1*2] - t[0+0*2] );
  }
/*
  Normalize the first three functions.
*/
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      phi[i+j*4]    = phi[i+j*4]    / area;
      dphidx[i+j*4] = dphidx[i+j*4] / area;
      dphidy[i+j*4] = dphidy[i+j*4] / area;
    }
  }
/*
  Compute the cubic bubble function.
*/
  for ( j = 0; j < n; j++ )
  {
       phi[3+j*4] = 27.0 * phi[0+j*4] * phi[1+j*4] * phi[2+j*4];

    dphidx[3+j*4] = 27.0 * (
                    dphidx[0+j*4] *    phi[1+j*4] *    phi[2+j*4]
                    +  phi[0+j*4] * dphidx[1+j*4] *    phi[2+j*4]
                    +  phi[0+j*4] *    phi[1+j*4] * dphidx[2+j*4] );

    dphidy[3+j*4] = 27.0 * (
                    dphidy[0+j*4] *    phi[1+j*4] *    phi[2+j*4]
                    +  phi[0+j*4] * dphidy[1+j*4] *    phi[2+j*4]
                    +  phi[0+j*4] *    phi[1+j*4] * dphidy[2+j*4] );
  }
/*
  Subtract 1/3 of the cubic bubble function from each of the three linears.
*/
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
         phi[i+j*4] =    phi[i+j*4] -    phi[3+j*4] / 3.0;
      dphidx[i+j*4] = dphidx[i+j*4] - dphidx[3+j*4] / 3.0;
      dphidy[i+j*4] = dphidy[i+j*4] - dphidy[3+j*4] / 3.0;
    }
  }

  return;
}
/******************************************************************************/

void basis_mn_t4_test ( )

/******************************************************************************/
/*
  Purpose:

    BASIS_MN_T4_TEST verifies BASIS_MN_T4.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 February 2006

  Author:

    John Burkardt

  Parameters:
*/
{
# define NODE_NUM 4

  double dphidx[NODE_NUM*NODE_NUM];
  double dphidy[NODE_NUM*NODE_NUM];
  int i;
  int j;
  double phi[NODE_NUM*NODE_NUM];
  double sum_x;
  double sum_y;
  double t[2*NODE_NUM] = { 
    2.0, 0.0,
    4.0, 2.0,
    0.0, 4.0,
    2.0, 2.0 };

  printf ( "\n" );
  printf ( "BASIS_MN_T4_TEST:\n" );
  printf ( "  Verify basis functions for element T4.\n" );
  printf ( "\n" );
  printf ( "  Number of nodes = %d\n", NODE_NUM );

  printf ( "\n" );
  printf ( "  Physical Nodes:\n" );
  printf ( "\n" );
  for ( j = 0; j < NODE_NUM; j++ )
  {
    printf ( "  %3d  %10g  %10g\n", j, t[0+j*2], t[1+j*2] );
  }

  printf ( "\n" );
  printf ( "  The basis function values at basis nodes\n" );
  printf ( "  should form the identity matrix.\n" );
  printf ( "\n" );

  basis_mn_t4 ( t, NODE_NUM, t, phi, dphidx, dphidy );

  for ( j = 0; j < NODE_NUM; j++ )
  {
    for ( i = 0; i < NODE_NUM; i++ )
    {
      printf ( "  %10g", phi[i+j*NODE_NUM] );
    }
    printf ( "\n" );
  }

  printf ( "\n" );
  printf ( "  The X and Y derivatives should sum to 0.\n" );
  printf ( "\n" );
  printf ( "  dPhidX sum, dPhidY sum:\n" );
  printf ( "\n" );

  for ( j = 0; j < NODE_NUM; j++ )
  {
    sum_x = 0.0;
    sum_y = 0.0;
    for ( i = 0; i < NODE_NUM; i++ )
    {
      sum_x = sum_x + dphidx[i+j*NODE_NUM];
      sum_y = sum_y + dphidy[i+j*NODE_NUM];
    }
    printf ( "  %10g  %10g\n", sum_x, sum_y );
  }
  return;
# undef NODE_NUM
}
/******************************************************************************/

void basis_mn_t6 ( double t[2*6], int n, double p[], double phi[], 
  double dphidx[], double dphidy[] )

/******************************************************************************/
/*
  Purpose:

    BASIS_MN_T6: all bases at N points for a T6 element.

  Discussion:

    The routine is given the coordinates of the vertices and midside
    nodes of a triangle.  It works directly with these coordinates, and does 
    not refer to a reference element.

    This routine requires that the midside nodes be "in line"
    with the vertices, that is, that the sides of the triangle be
    straight.  However, the midside nodes do not actually have to
    be halfway along the side of the triangle.  

  The physical element T6:

    This picture indicates the assumed ordering of the six nodes
    of the triangle.

    |
    |   
    |        3
    |       / \
    |      /   \
    Y     6     5
    |    /       \
    |   /         \
    |  1-----4-----2
    |
    +--------X-------->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 February 2006

  Author:

    John Burkardt

  Parameters:

    Input, double T[2*6], the nodal oordinates of the element.
    It is common to list these points in counter clockwise order.

    Input, int N, the number of evaluation points.

    Input, double P[2*N], the coordinates of the point where
    the basis functions are to be evaluated.

    Output, double PHI[6*N], the value of the basis functions at P.

    Output, double DPHIDX[6*N], DPHIDY[6*N], the value of the X 
    and Y derivatives of the basis functions at P.
*/
{
  double gn;
  double gx;
  double hn;
  double hx;
  int j;

  for ( j = 0; j < n; j++ )
  {
/*
  Basis function 1: PHI(X,Y) = G(3,2) * H(6,4) / normalization.
*/
    gx = ( p[0+j*2] - t[0+1*2] ) * ( t[1+2*2] - t[1+1*2] ) 
       - ( t[0+2*2] - t[0+1*2] ) * ( p[1+j*2] - t[1+1*2] );

    gn = ( t[0+0*2] - t[0+1*2] ) * ( t[1+2*2] - t[1+1*2] ) 
       - ( t[0+2*2] - t[0+1*2] ) * ( t[1+0*2] - t[1+1*2] );

    hx = ( p[0+j*2] - t[0+3*2] ) * ( t[1+5*2] - t[1+3*2] ) 
       - ( t[0+5*2] - t[0+3*2] ) * ( p[1+j*2] - t[1+3*2] );

    hn = ( t[0+0*2] - t[0+3*2] ) * ( t[1+5*2] - t[1+3*2] ) 
       - ( t[0+5*2] - t[0+3*2] ) * ( t[1+0*2] - t[1+3*2] );

    phi[0+j*6] =     ( gx * hx ) / ( gn * hn );
    dphidx[0+j*6] =  (      ( t[1+2*2] - t[1+1*2] ) * hx 
                     + gx * ( t[1+5*2] - t[1+3*2] ) ) / ( gn * hn );
    dphidy[0+j*6] = -(      ( t[0+2*2] - t[0+1*2] ) * hx 
                     + gx * ( t[0+5*2] - t[0+3*2] ) ) / ( gn * hn );
/*
  Basis function 2: PHI(X,Y) = G(3,1) * H(4,5) / normalization.
*/
    gx = ( p[0+j*2] - t[0+0*2] ) * ( t[1+2*2] - t[1+0*2] ) 
       - ( t[0+2*2] - t[0+0*2] ) * ( p[1+j*2] - t[1+0*2] );

    gn = ( t[0+1*2] - t[0+0*2] ) * ( t[1+2*2] - t[1+0*2] ) 
       - ( t[0+2*2] - t[0+0*2] ) * ( t[1+1*2] - t[1+0*2] );

    hx = ( p[0+j*2] - t[0+4*2] ) * ( t[1+3*2] - t[1+4*2] ) 
       - ( t[0+3*2] - t[0+4*2] ) * ( p[1+j*2] - t[1+4*2] );

    hn = ( t[0+1*2] - t[0+4*2] ) * ( t[1+3*2] - t[1+4*2] ) 
       - ( t[0+3*2] - t[0+4*2] ) * ( t[1+1*2] - t[1+4*2] );

    phi[1+j*6] = ( gx * hx ) / ( gn * hn );
    dphidx[1+j*6] =  (      ( t[1+2*2] - t[1+0*2] ) * hx 
                     + gx * ( t[1+3*2] - t[1+4*2] ) ) / ( gn * hn );
    dphidy[1+j*6] = -(      ( t[0+2*2] - t[0+0*2] ) * hx 
                     + gx * ( t[0+3*2] - t[0+4*2] ) ) / ( gn * hn );
/*
  Basis function 3: PHI(X,Y) = G(1,2) * H(5,6) / normalization.
*/
    gx = ( p[0+j*2] - t[0+1*2] ) * ( t[1+0*2] - t[1+1*2] ) 
       - ( t[0+0*2] - t[0+1*2] ) * ( p[1+j*2] - t[1+1*2] );

    gn = ( t[0+2*2] - t[0+1*2] ) * ( t[1+0*2] - t[1+1*2] ) 
       - ( t[0+0*2] - t[0+1*2] ) * ( t[1+2*2] - t[1+1*2] );

    hx = ( p[0+j*2] - t[0+5*2] ) * ( t[1+4*2] - t[1+5*2] ) 
       - ( t[0+4*2] - t[0+5*2] ) * ( p[1+j*2] - t[1+5*2] );

    hn = ( t[0+2*2] - t[0+5*2] ) * ( t[1+4*2] - t[1+5*2] ) 
       - ( t[0+4*2] - t[0+5*2] ) * ( t[1+2*2] - t[1+5*2] );

    phi[2+j*6] = ( gx * hx ) / ( gn * hn );
    dphidx[2+j*6] =  (      ( t[1+0*2] - t[1+1*2] ) * hx 
                     + gx * ( t[1+4*2] - t[1+5*2] ) ) / ( gn * hn );
    dphidy[2+j*6] = -(      ( t[0+0*2] - t[0+1*2] ) * hx 
                     + gx * ( t[0+4*2] - t[0+5*2] ) ) / ( gn * hn );
/*
  Basis function 4: PHI(X,Y) = G(1,3) * H(2,3) / normalization.
*/
    gx = ( p[0+j*2] - t[0+2*2] ) * ( t[1+0*2] - t[1+2*2] ) 
       - ( t[0+0*2] - t[0+2*2] ) * ( p[1+j*2] - t[1+2*2] );

    gn = ( t[0+3*2] - t[0+2*2] ) * ( t[1+0*2] - t[1+2*2] ) 
       - ( t[0+0*2] - t[0+2*2] ) * ( t[1+3*2] - t[1+2*2] );

    hx = ( p[0+j*2] - t[0+2*2] ) * ( t[1+1*2] - t[1+2*2] ) 
       - ( t[0+1*2] - t[0+2*2] ) * ( p[1+j*2] - t[1+2*2] );

    hn = ( t[0+3*2] - t[0+2*2] ) * ( t[1+1*2] - t[1+2*2] ) 
       - ( t[0+1*2] - t[0+2*2] ) * ( t[1+3*2] - t[1+2*2] );

    phi[3+j*6] = ( gx * hx ) / ( gn * hn );
    dphidx[3+j*6] =  (      ( t[1+0*2] - t[1+2*2] ) * hx 
                     + gx * ( t[1+1*2] - t[1+2*2] ) ) / ( gn * hn );
    dphidy[3+j*6] = -(      ( t[0+0*2] - t[0+2*2] ) * hx 
                     + gx * ( t[0+1*2] - t[0+2*2] ) ) / ( gn * hn );
/*
  Basis function 5: PHI(X,Y) = G(2,1) * H(3,1) / normalization.
*/
    gx = ( p[0+j*2] - t[0+0*2] ) * ( t[1+1*2] - t[1+0*2] ) 
       - ( t[0+1*2] - t[0+0*2] ) * ( p[1+j*2] - t[1+0*2] );

    gn = ( t[0+4*2] - t[0+0*2] ) * ( t[1+1*2] - t[1+0*2] ) 
       - ( t[0+1*2] - t[0+0*2] ) * ( t[1+4*2] - t[1+0*2] );

    hx = ( p[0+j*2] - t[0+0*2] ) * ( t[1+2*2] - t[1+0*2] ) 
       - ( t[0+2*2] - t[0+0*2] ) * ( p[1+j*2] - t[1+0*2] );

    hn = ( t[0+4*2] - t[0+0*2] ) * ( t[1+2*2] - t[1+0*2] ) 
       - ( t[0+2*2] - t[0+0*2] ) * ( t[1+4*2] - t[1+0*2] );

    phi[4+j*6] = ( gx * hx ) / ( gn * hn );
    dphidx[4+j*6] =  (      ( t[1+1*2] - t[1+0*2] ) * hx 
                     + gx * ( t[1+2*2] - t[1+0*2] ) ) / ( gn * hn );
    dphidy[4+j*6] = -(      ( t[0+1*2] - t[0+0*2] ) * hx 
                     + gx * ( t[0+2*2] - t[0+0*2] ) ) / ( gn * hn );
/*
  Basis function 6: PHI(X,Y) = G(1,2) * H(3,2) / normalization.
*/
    gx = ( p[0+j*2] - t[0+1*2] ) * ( t[1+0*2] - t[1+1*2] ) 
       - ( t[0+0*2] - t[0+1*2] ) * ( p[1+j*2] - t[1+1*2] );

    gn = ( t[0+5*2] - t[0+1*2] ) * ( t[1+0*2] - t[1+1*2] ) 
       - ( t[0+0*2] - t[0+1*2] ) * ( t[1+5*2] - t[1+1*2] );

    hx = ( p[0+j*2] - t[0+1*2] ) * ( t[1+2*2] - t[1+1*2] ) 
       - ( t[0+2*2] - t[0+1*2] ) * ( p[1+j*2] - t[1+1*2] );

    hn = ( t[0+5*2] - t[0+1*2] ) * ( t[1+2*2] - t[1+1*2] ) 
       - ( t[0+2*2] - t[0+1*2] ) * ( t[1+5*2] - t[1+1*2] );

    phi[5+j*6] = ( gx * hx ) / ( gn * hn );
    dphidx[5+j*6] =  (      ( t[1+0*2] - t[1+1*2] ) * hx 
                     + gx * ( t[1+2*2] - t[1+1*2] ) ) / ( gn * hn );
    dphidy[5+j*6] = -(      ( t[0+0*2] - t[0+1*2] ) * hx 
                     + gx * ( t[0+2*2] - t[0+1*2] ) ) / ( gn * hn );
  }

  return;
}
/******************************************************************************/

void basis_mn_t6_test ( )

/******************************************************************************/
/*
  Purpose:

    BASIS_MN_T6_TEST verifies BASIS_MN_T6.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 February 2006

  Author:

    John Burkardt

  Parameters:

    None
*/
{
# define NODE_NUM 6

  double dphidx[NODE_NUM*NODE_NUM];
  double dphidy[NODE_NUM*NODE_NUM];
  int i;
  int j;
  double phi[NODE_NUM*NODE_NUM];
  double sum_x;
  double sum_y;
  double t[2*NODE_NUM] = {
    2.0, 0.0,
    4.0, 3.0,
    0.0, 4.0,
    3.0, 1.5,
    2.0, 3.5,
    1.0, 2.0 };

  printf ( "\n" );
  printf ( "BASIS_MN_T6_TEST:\n" );
  printf ( "  Verify basis functions for element T6.\n" );
  printf ( "\n" );
  printf ( "  Number of nodes = %d\n", NODE_NUM );

  printf ( "\n" );
  printf ( "  Physical Nodes:\n" );
  printf ( "\n" );
  for ( j = 0; j < NODE_NUM; j++ )
  {
    printf ( "  %3d  %10g  %10g\n", j, t[0+j*2], t[1+j*2] );
  }
 
  printf ( "\n" );
  printf ( "  The basis function values at basis nodes\n" );
  printf ( "  should form the identity matrix.\n" );
  printf ( "\n" );

  basis_mn_t6 ( t, NODE_NUM, t, phi, dphidx, dphidy );

  for ( i = 0; i < NODE_NUM; i++ )
  {
    for ( j = 0; j < NODE_NUM; j++ )
    {
      printf ( "  %7g", phi[i+j*NODE_NUM] );
    }
    printf ( "\n" );
  }

  printf ( "\n" );
  printf ( "  The X and Y derivatives should sum to 0.\n" );
  printf ( "\n" );
  printf ( "  dPhidX sum, dPhidY sum:\n" );
  printf ( "\n" );
  for ( j = 0; j < NODE_NUM; j++ )
  {
    sum_x = 0.0;
    for ( i = 0; i < NODE_NUM; i++ )
    {
      sum_x = sum_x + dphidx[i+j*NODE_NUM];
    }
    sum_y = 0.0;
    for ( i = 0; i < NODE_NUM; i++ )
    {
      sum_y = sum_y + dphidy[i+j*NODE_NUM];
    }
    printf ( "  %10g  %10g\n", sum_x, sum_y );
  }

  return;
# undef NODE_NUM
}
/******************************************************************************/

char ch_cap ( char ch )

/******************************************************************************/
/*
  Purpose:

    CH_CAP capitalizes a single character.

  Discussion:

    This routine should be equivalent to the library "toupper" function.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 July 1998

  Author:

    John Burkardt

  Parameters:

    Input, char CH, the character to capitalize.

    Output, char CH_CAP, the capitalized character.
*/
{
  if ( 97 <= ch && ch <= 122 )
  {
    ch = ch - 32;
  }

  return ch;
}
/******************************************************************************/

double degrees_to_radians ( double angle )

/******************************************************************************/
/*
  Purpose: 

    DEGREES_TO_RADIANS converts an angle from degrees to radians.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, double ANGLE, an angle in degrees.

    Output, double DEGREES_TO_RADIANS, the equivalent angle
    in radians.
*/
{
# define PI 3.141592653589793

  return ( angle * PI / 180.0 );

# undef PI
}
/******************************************************************************/

void derivative_average_t3 ( int node_num, double node_xy[], int element_num,
  int element_node[], double c[], double dcdx[], double dcdy[] )

/******************************************************************************/
/*
  Purpose:

    DERIVATIVE_AVERAGE_T3 averages derivatives at the nodes of a T3 mesh.

  Discussion:

    This routine can be used to compute an averaged nodal value of any
    quantity associated with the finite element function.  At a node
    that is shared by several elements, the fundamental function
    U will be continuous, but its spatial derivatives, for instance,
    will generally be discontinuous.  This routine computes the
    value of the spatial derivatives in each element, and averages
    them, to make a reasonable assignment of a nodal value.

    Note that the ELEMENT_NODE array is assumed to be 1-based, rather
    than 0-based.  Thus, entries from this array must be decreased by
    1 before being used!

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NODE_NUM, the number of nodes.

    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.

    Input, int ELEMENT_NUM, the number of elements.

    Input, int ELEMENT_NODE[3*ELEMENT_NUM], the element->node data.

    Input, double C[NODE_NUM], the finite element coefficient vector.

    Output, double DCDX[NODE_NUM], DCDY[NODE_NUM], the averaged
    values of dCdX and dCdY at the nodes.
*/
{
# define OFFSET 1

  int dim;
  double dphidx[3*3];
  double dphidy[3*3];
  int element;
  int j;
  int node;
  int node_count[node_num];
  int node_global1;
  int node_global2;
  int node_local1;
  int node_local2;
  double phi[3*3];
  double t[2*3];

  for ( node = 0; node < node_num; node++ )
  {
    node_count[node] = 0;
    dcdx[node] = 0.0;
    dcdy[node] = 0.0;
  }
/*
  Consider every element.
*/
  for ( element = 0; element < element_num; element++ )
  {
/*
  Get the coordinates of the nodes of the element.
*/
    for ( j = 0; j < 3; j++ )
    {
      for ( dim = 0; dim < 2; dim++ )
      {
        t[dim+2*j] = node_xy[dim+(element_node[j+element*3]-OFFSET)];
      }
    }
/*
  Evaluate the X and Y derivatives of the 3 basis functions at the
  3 nodes.
*/
    basis_mn_t3 ( t, 3, t, phi, dphidx, dphidy );
/*
  Evaluate dCdX and dCdY at each node in the element, and add
  them to the running totals.
*/
    for ( node_local1 = 0; node_local1 < 3; node_local1++ )
    {
      node_global1 = element_node[node_local1+element*3]-OFFSET;

      for ( node_local2 = 0; node_local2 < 3; node_local2++ )
      {
        node_global2 = element_node[node_local2+element*3]-OFFSET;

        dcdx[node_global1] = dcdx[node_global1]
          + c[node_global2] * dphidx[node_local2+node_local1*3];

        dcdy[node_global1] = dcdy[node_global1] 
          + c[node_global2] * dphidy[node_local2+node_local1*3];
      }
      node_count[node_global1] = node_count[node_global1] + 1;
    }
  }
/*
  Average the running totals.
*/
  for ( node = 0; node < node_num; node++ )
  {
    dcdx[node] = dcdx[node] / ( double ) node_count[node];
    dcdy[node] = dcdy[node] / ( double ) node_count[node];
  }
  return;
# undef OFFSET
}
/******************************************************************************/

void div_q4 ( int m, int n, double u[], double v[], double xlo, double xhi, 
  double ylo, double yhi, double div[], double vort[] )

/******************************************************************************/
/*
  Purpose: 

    DIV_Q4 estimates the divergence and vorticity of a discrete field.

  Discussion:

    The routine is given the values of a vector field ( U(X,Y), V(X,Y) ) at
    an array of points ( X(1:M), Y(1:N) ).

    The routine models the vector field over the interior of this region using
    a bilinear interpolant.  It then uses the interpolant to estimate the
    value of the divergence:

      DIV(X,Y) = dU/dX + dV/dY

    and the vorticity:

      VORT(X,Y) = dV/dX - dU/dY

    at the center point of each of the bilinear elements.

        |       |       |
      (3,1)---(3,2)---(3,3)---
        |       |       |
        | [2,1] | [2,2] |
        |       |       |
      (2,1)---(2,2)---(2,3)---
        |       |       |
        | [1,1] | [1,2] |
        |       |       |
      (1,1)---(1,2)---(1,3)---

    Here, the nodes labeled with parentheses represent the points at
    which the original (U,V) data is given, while the nodes labeled
    with square brackets represent the centers of the bilinear
    elements, where the approximations to the divergence and vorticity
    are made.

    The reason for evaluating the divergence and vorticity in this way
    is that the bilinear interpolant to the (U,V) data is not
    differentiable at the boundaries of the elements, nor especially at
    the nodes, but is an (infinitely differentiable) bilinear function
    in the interior of each element.  If a value at the original nodes
    is strongly desired, then the average at the four surrounding
    central nodes may be taken.

  Element Q4:

    |
    1  4-----3
    |  |     |
    |  |     |
    S  |     |
    |  |     |
    |  |     |
    0  1-----2
    |
    +--0--R--1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of data rows.  M must be at least 2.

    Input, int N, the number of data columns.  N must be at least 2.

    Input, double U[M*N], V[M*N], the value of the components 
    of a vector quantity whose divergence and vorticity are desired. 
    A common example would be that U and V are the horizontal and 
    vertical velocity component of a flow field.

    Input, double XLO, XHI, the minimum and maximum X coordinates.

    Input, double YLO, YHI, the minimum and maximum Y coordinates.

    Output, double DIV[(M-1)*(N-1)], an estimate for
    the divergence in the bilinear element that lies between
    data rows I and I+1, and data columns J and J+1.

    Output, double VORT[(M-1)*(N-1)], an estimate for
    the vorticity in the bilinear element that lies between
    data rows I and I+1, and data columns J and J+1.
*/
{
  double dphidx[4];
  double dphidy[4];
  int i;
  int j;
  double p[2];
  double phi[4];
  double q[2*4];
  double xl;
  double xr;
  double yb;
  double yt;

  if ( m <= 1 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "DIV_Q4 - Fatal error!\n" );
    fprintf ( stderr, "  M must be at least 2,\n" );
    fprintf ( stderr, "  but the input value of M is %d\n", m );
    exit ( 1 );
  }

  if ( n <= 1 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "DIV_Q4 - Fatal error!\n" );
    fprintf ( stderr, "  N must be at least 2,\n" );
    fprintf ( stderr, "  but the input value of N is %d\n", n );
    exit ( 1 );
  }

  if ( xhi == xlo )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "DIV_Q4 - Fatal error!\n" );
    fprintf ( stderr, "  XHI and XLO must be distinct,\n" );
    fprintf ( stderr, "  but the input value of XLO is %g\n", xlo );
    fprintf ( stderr, "  and the input value of XHI is %g\n", xhi );
    exit ( 1 );
  }

  if ( yhi == ylo )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "DIV_Q4 - Fatal error!\n" );
    fprintf ( stderr, "  YHI and YLO must be distinct,\n" );
    fprintf ( stderr, "  but the input value of YLO is %g\n", ylo );
    fprintf ( stderr, "  and the input value of YHI is %g\n", yhi );
    exit ( 1 );
  }

  for ( i = 1; i <= m-1; i++ )
  {
    yb = ( ( double ) ( 2 * m - 2 * i     ) * ylo   
         + ( double ) (         2 * i - 2 ) * yhi ) 
         / ( double ) ( 2 * m         - 2 );
    p[1] = ( ( double ) ( 2 * m - 2 * i - 1 ) * ylo   
           + ( double ) (         2 * i - 1 ) * yhi ) 
           / ( double ) ( 2 * m         - 2 );
    yt = ( ( double ) ( 2 * m - 2 * i - 2 ) * ylo   
         + ( double ) (         2 * i     ) * yhi ) 
         / ( double ) ( 2 * m         - 2 );

    q[1+0*2] = yb;
    q[1+1*2] = yb;
    q[1+2*2] = yt;
    q[1+3*2] = yt;

    for ( j = 1; j <= n-1; j++ )
    {
      xl = ( ( double ) ( 2 * n - 2 * j     ) * xlo   
           + ( double ) (         2 * j - 2 ) * xhi ) 
           / ( double ) ( 2 * n         - 2 );
      p[0] = ( ( double ) ( 2 * n - 2 * j - 1 ) * xlo   
             + ( double ) (         2 * j - 1 ) * xhi ) 
             / ( double ) ( 2 * n         - 2 );
      xr = ( ( double ) ( 2 * n - 2 * j - 2 ) * xlo   
           + ( double ) (         2 * j     ) * xhi ) 
           / ( double ) ( 2 * n         - 2 );

      q[0+0*2] = xl;
      q[0+1*2] = xr;
      q[0+2*2] = xr;
      q[0+3*2] = xl;
/*
  Evaluate the basis function and derivatives at the center of the element.
*/
      basis_mn_q4 ( q, 1, p, phi, dphidx, dphidy );
/*
  Note the following formula for the value of U and V at the same
  point that the divergence and vorticity are being evaluated.

         umid =  u(i  ,j  ) * phi[0] &
               + u(i  ,j+1) * phi[1] &
               + u(i+1,j+1) * phi[2] &
               + u(i+1,j  ) * phi[3] 

         vmid =  v(i  ,j  ) * phi[0] &
               + v(i  ,j+1) * phi[1] &
               + v(i+1,j+1) * phi[2] &
               + v(i+1,j  ) * phi[3] 
*/
      div[i-1+(j-1)*(m-1)] =  
                  u[i-1+(j-1)*m] * dphidx[0] + v[i-1+(j-1)*m] * dphidy[0] 
                + u[i-1+(j  )*m] * dphidx[1] + v[i-1+(j  )*m] * dphidy[1] 
                + u[i  +(j  )*m] * dphidx[2] + v[i  +(j  )*m] * dphidy[2] 
                + u[i  +(j-1)*m] * dphidx[3] + v[i  +(j-1)*m] * dphidy[3];

      vort[i-1+(j-1)*(m-1)] =  
                   v[i-1+(j-1)*m] * dphidx[0] - u[i-1+(j-1)*m] * dphidy[0] 
                 + v[i-1+(j  )*m] * dphidx[1] - u[i-1+(j  )*m] * dphidy[1] 
                 + v[i  +(j  )*m] * dphidx[2] - u[i  +(j  )*m] * dphidy[2] 
                 + v[i  +(j-1)*m] * dphidx[3] - u[i  +(j-1)*m] * dphidy[3]; 
    }
  }

  return;
}
/******************************************************************************/

char *element_code ( int i )

/******************************************************************************/
/*
  Purpose:

    ELEMENT_CODE returns the code for each element.

  List:

    I  ELEMENT_CODE   Definition
    -  ------------   ----------
    1  Q4             4 node linear Lagrange/serendipity quadrilateral;
    2  Q8             8 node quadratic serendipity quadrilateral;
    3  Q9             9 node quadratic Lagrange quadrilateral;
    4  Q12            12 node cubic serendipity quadrilateral;
    5  Q16            16 node cubic Lagrange quadrilateral;
    6  QL             6 node linear/quadratic quadrilateral;
    7  T3             3 node linear triangle;
    8  T4             4 node cubic bubble triangle
    9  T6             6 node quadratic triangle;
   10  T10            10 node cubic triangle.
 
  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int I, the number of the element.  

    Output, char *ELEMENT_CODE, the code for the element.
*/
{
  char *value;

  if ( i == 1 )
  {
    value = ( char * ) malloc ( 3 * sizeof ( char ) );
    strcpy ( value, "Q4" );
  }
  else if ( i == 2 )
  {
    value = ( char * ) malloc ( 3 * sizeof ( char ) );
    strcpy ( value, "Q8" );
  }
  else if ( i == 3 )
  {
    value = ( char * ) malloc ( 3 * sizeof ( char ) );
    strcpy ( value, "Q9" );
  }
  else if ( i == 4 )
  {
    value = ( char * ) malloc ( 4 * sizeof ( char ) );
    strcpy ( value, "Q12" );
  }
  else if ( i == 5 )
  {
    value = ( char * ) malloc ( 4 * sizeof ( char ) );
    strcpy ( value, "Q16" );
  }
  else if ( i == 6 )
  {
    value = ( char * ) malloc ( 3 * sizeof ( char ) );
    strcpy ( value, "QL" );
  }
  else if ( i == 7 )
  {
    value = ( char * ) malloc ( 3 * sizeof ( char ) );
    strcpy ( value, "T3" );
  }
  else if ( i == 8 )
  {
    value = ( char * ) malloc ( 3 * sizeof ( char ) );
    strcpy ( value, "T4" );
  }
  else if ( i == 9 )
  {
    value = ( char * ) malloc ( 3 * sizeof ( char ) );
    strcpy ( value, "T6" );
  }
  else if ( i == 10 )
  {
    value = ( char * ) malloc ( 4 * sizeof ( char ) );
    strcpy ( value, "T10" );
  }
  else
  {
    value = ( char * ) malloc ( 4 * sizeof ( char ) );
    strcpy ( value, "???" );
  }

  return value;
}
/******************************************************************************/

void elements_eps ( char *file_name, int node_num, double node_xy[], char *code, 
  int element_num, int element_mask[], int element_node[], int node_show, 
  int element_show )

/******************************************************************************/
/*
  Purpose: 

    ELEMENTS_EPS creates an EPS file image of the elements of a grid.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, char *FILE_NAME, the name of the file to create.

    Input, int NODE_NUM, the number of nodes.

    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.

    Input, char *CODE, the code for the element.

    Input, int ELEMENT_NUM, the number of elements.

    Input, int ELEMENT_MASK[ELEMENT_NUM], a mask for the elements.
    Only elements with a TRUE mask will be shown.

    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], the nodes making up
    each element.

    Input, int NODE_SHOW:
    0, do not show nodes;
    1, show nodes;
    2, show nodes and label them.

    Input, int TRIANGLE_SHOW:
    0, do not show triangles;
    1, show triangles;
    2, show triangles and label them.
*/
{
  double ave_x;
  double ave_y;
  int circle_size = 3;
  int delta;
  int element;
  int element_order;
  FILE *file_unit;
  int local;
  int node;
  int *node_mask;
  double x_max;
  double x_min;
  int x_ps;
  int x_ps_max = 576;
  int x_ps_max_clip = 594;
  int x_ps_min = 36;
  int x_ps_min_clip = 18;
  double x_scale;
  double y_max;
  double y_min;
  int y_ps;
  int y_ps_max = 666;
  int y_ps_max_clip = 684;
  int y_ps_min = 126;
  int y_ps_min_clip = 108;
  double y_scale;

  element_order = order_code ( code );
/*
  Determine which nodes are visible, controlled by which elements are visible.
*/
  node_mask = ( int * ) malloc ( node_num * sizeof ( int ) );
  for ( node = 0; node < node_num; node++ )
  {
    node_mask[node] = 0;
  }
  for ( element = 0; element < element_num; element++ )
  {
    if ( element_mask[element] )
    {
      for ( local = 0; local < element_order; local++ )
      {
        node = element_node[local+element*element_order]-1;
        node_mask[node] = 1;
      }
    }
  }
/*
  We need to do some figuring here, so that we can determine
  the range of the data, and hence the height and width
  of the piece of paper.
*/
  x_max = -r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
     if ( node_mask[node] )
     {
       if ( x_max < node_xy[0+node*2] )
       {
         x_max = node_xy[0+node*2];
       }
    }
  }
  x_min = r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
    if ( node_mask[node] )
    {
       if ( node_xy[0+node*2] < x_min )
       {
         x_min = node_xy[0+node*2];
       }
    }
  }
  x_scale = x_max - x_min;

  x_max = x_max + 0.05 * x_scale;
  x_min = x_min - 0.05 * x_scale;
  x_scale = x_max - x_min;

  y_max = -r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
    if ( node_mask[node] )
    {
      if ( y_max < node_xy[1+node*2] )
      {
        y_max = node_xy[1+node*2];
      }
    }
  }
  y_min = r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
    if ( node_mask[node] )
    {
      if ( node_xy[1+node*2] < y_min )
      {
        y_min = node_xy[1+node*2];
      }
    }
  }
  y_scale = y_max - y_min;

  y_max = y_max + 0.05 * y_scale;
  y_min = y_min - 0.05 * y_scale;
  y_scale = y_max - y_min;

  if ( x_scale < y_scale )
  {
    delta = r8_nint ( ( double ) ( x_ps_max - x_ps_min )
      * ( y_scale - x_scale ) / ( 2.0 * y_scale ) );

    x_ps_max = x_ps_max - delta;
    x_ps_min = x_ps_min + delta;

    x_ps_max_clip = x_ps_max_clip - delta;
    x_ps_min_clip = x_ps_min_clip + delta;

    x_scale = y_scale;
  }
  else if ( y_scale < x_scale )
  {
    delta = r8_nint ( ( double ) ( y_ps_max - y_ps_min ) 
      * ( x_scale - y_scale ) / ( 2.0 * x_scale ) );

    y_ps_max = y_ps_max - delta;
    y_ps_min = y_ps_min + delta;

    y_ps_max_clip = y_ps_max_clip - delta;
    y_ps_min_clip = y_ps_min_clip + delta;

    y_scale = x_scale;
  }

  file_unit = fopen ( file_name, "wt" );

  if ( !file_unit )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "ELEMENTS_EPS - Fatal error!\n" );
    fprintf ( stderr, "  Could not open the output EPS file.\n" );
    exit ( 1 );
  }

  fprintf ( file_unit, "%%!PS-Adobe-3.0 EPSF-3.0\n" );
  fprintf ( file_unit, "%%Creator: elements_eps.C\n" );
  fprintf ( file_unit, "%%Title: %s\n", file_name );

  fprintf ( file_unit, "%%Pages: 1\n" );
  fprintf ( file_unit, "%%BoundingBox:  %d  %d  %d  %d\n", 
    x_ps_min, y_ps_min, x_ps_max, y_ps_max );
  fprintf ( file_unit, "%%Document-Fonts: Times-Roman\n" );
  fprintf ( file_unit, "%%LanguageLevel: 1\n" );
  fprintf ( file_unit, "%%EndComments\n" );
  fprintf ( file_unit, "%%BeginProlog\n" );
  fprintf ( file_unit, "/inch {72 mul} def\n" );
  fprintf ( file_unit, "%%EndProlog\n" );
  fprintf ( file_unit, "%%Page:      1     1\n" );
  fprintf ( file_unit, "save\n" );
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, "%% Set the RGB line color to very light gray.\n" );
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, " 0.9000 0.9000 0.9000 setrgbcolor\n" );
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, "%% Draw a gray border around the page.\n" );
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, "newpath\n" );
  fprintf ( file_unit, "%d  %d  moveto\n", x_ps_min, y_ps_min );
  fprintf ( file_unit, "%d  %d  lineto\n", x_ps_max, y_ps_min );
  fprintf ( file_unit, "%d  %d  lineto\n", x_ps_max, y_ps_max );
  fprintf ( file_unit, "%d  %d  lineto\n", x_ps_min, y_ps_max );
  fprintf ( file_unit, "%d  %d  lineto\n", x_ps_min, y_ps_min );
  fprintf ( file_unit, "stroke\n" );
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, "%% Set RGB line color to black.\n" );
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, " 0.0000 0.0000 0.0000 setrgbcolor\n" );
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, "%%  Set the font and its size:\n" );
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, "/Times-Roman findfont\n" );
  fprintf ( file_unit, "0.50 inch scalefont\n" );
  fprintf ( file_unit, "setfont\n" );
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, "%%  Print a title:\n" );
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, "%%  210  702 moveto\n" );
  fprintf ( file_unit, "%%(Pointset) show\n" );
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, "%% Define a clipping polygon\n" );
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, "newpath\n" );
  fprintf ( file_unit, "%d  %d  moveto\n", x_ps_min_clip, y_ps_min_clip );
  fprintf ( file_unit, "%d  %d  lineto\n", x_ps_max_clip, y_ps_min_clip );
  fprintf ( file_unit, "%d  %d  lineto\n", x_ps_max_clip, y_ps_max_clip );
  fprintf ( file_unit, "%d  %d  lineto\n", x_ps_min_clip, y_ps_max_clip );
  fprintf ( file_unit, "%d  %d  lineto\n", x_ps_min_clip, y_ps_min_clip );
  fprintf ( file_unit, "clip newpath\n" );
/*
  Draw the nodes.
*/
  if ( 1 <= node_show )
  {
    fprintf ( file_unit, "%%\n" );
    fprintf ( file_unit, "%%  Draw filled dots at each node:\n" );
    fprintf ( file_unit, "%%\n" );
    fprintf ( file_unit, "%%  Set the color to blue:\n" );
    fprintf ( file_unit, "%%\n" );
    fprintf ( file_unit, "0.000  0.150  0.750  setrgbcolor\n" );
    fprintf ( file_unit, "%%\n" );

    for ( node = 0; node < node_num; node++ )
    {
      if ( node_mask[node] )
      {
        x_ps = ( int ) (
          ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )  
          + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) ) 
          / ( x_max                     - x_min ) );

        y_ps = ( int ) (
          ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )  
          + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) ) 
          / ( y_max                     - y_min ) );

        fprintf ( file_unit, "newpath  %d  %d  %d  0 360 arc closepath fill\n",
          x_ps, y_ps, circle_size );
      }
    }
  }
/*
  Label the nodes.
*/
  if ( 2 <= node_show )
  {
    fprintf ( file_unit, "%%\n" );
    fprintf ( file_unit, "%%  Label the nodes:\n" );
    fprintf ( file_unit, "%%\n" );
    fprintf ( file_unit, "%%  Set the color to darker blue:\n" );
    fprintf ( file_unit, "%%\n" );
    fprintf ( file_unit, "0.000  0.250  0.850  setrgbcolor\n" );
    fprintf ( file_unit, "/Times-Roman findfont\n" );
    fprintf ( file_unit, "0.20 inch scalefont\n" );
    fprintf ( file_unit, "setfont\n" );

    fprintf ( file_unit, "%%\n" );

    for ( node = 0; node < node_num; node++ )
    { 
      if ( node_mask[node] )
      {
        x_ps = ( int ) (
          ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )  
          + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) ) 
          / ( x_max                     - x_min ) );

        y_ps = ( int ) (
          ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )  
          + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) ) 
          / ( y_max                     - y_min ) );

        fprintf ( file_unit, "newpath  %d  %d  moveto (%d) show\n",
          x_ps, y_ps + 5, node+1 );
      }
    }
  }
/*
  Draw the elements.
*/
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, "%%  Draw the element sides:\n" );
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, " 9.0000 0.0000 0.0000 setrgbcolor\n" );

  for ( element = 0; element < element_num; element++ )
  {
    if ( element_mask[element] )
    {
      local = 1;
      node = element_node[local-1+element*element_order] - 1;

      x_ps = ( int ) (
        ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )  
        + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) ) 
        / ( x_max                     - x_min ) );

      y_ps = ( int ) (
        ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )  
        + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) ) 
        / ( y_max                     - y_min ) );

      fprintf ( file_unit, "newpath  %d  %d\n", x_ps, y_ps );

      for ( ; ; )
      {
        local = next_boundary_node ( local, code );
        node = element_node[local-1+element*element_order] - 1;

        x_ps = ( int ) (
          ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )  
          + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) ) 
          / ( x_max                     - x_min ) );

        y_ps = ( int ) (
          ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )  
          + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) ) 
          / ( y_max                     - y_min ) );

        fprintf ( file_unit, "  %d  %d  lineto", x_ps, y_ps ); 

        if ( local == 1 )
        {
          break;
        }
      }
      fprintf ( file_unit, "stroke\n" );
    }
  }
/*
  Label the elements.
*/
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, "%%  Label the elements:\n" );
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, " 1.0000 0.0000 0.0000 setrgbcolor\n" );
  fprintf ( file_unit, "/Times-Roman findfont\n" );
  fprintf ( file_unit, "0.30 inch scalefont setfont\n" );

  for ( element = 0; element < element_num; element++ )
  {
    if ( element_mask[element] )
    {
      ave_x = 0.0;
      ave_y = 0.0;

      for ( local = 0; local < element_order; local++ )
      {
        node = element_node[local+element_order*element] - 1;

        ave_x = ave_x + node_xy[0+node*2];
        ave_y = ave_y + node_xy[1+node*2];
      }
      ave_x = ave_x / ( double ) ( element_order );
      ave_y = ave_y / ( double ) ( element_order );

      x_ps = ( int ) (
        ( ( x_max - ave_x         ) * ( double ) ( x_ps_min )  
        + (       + ave_x - x_min ) * ( double ) ( x_ps_max ) ) 
        / ( x_max         - x_min ) );

      y_ps = ( int ) (
        ( ( y_max - ave_y         ) * ( double ) ( y_ps_min )  
        + (         ave_y - y_min ) * ( double ) ( y_ps_max ) ) 
        / ( y_max         - y_min ) );

      fprintf ( file_unit, "newpath  %d  %d  moveto (%d) show\n", 
        x_ps, y_ps + 5, element+1 );
    }
  }
/*
  Finish up.
*/
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, "restore  showpage\n" );
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, "%%  End of page.\n" );
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, "%%Trailer\n" );
  fprintf ( file_unit, "%%EOF\n" );

  fclose ( file_unit );

  free ( node_mask );

  return;
}
/******************************************************************************/

int *grid_element ( char *code, int element_order, int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_ELEMENT returns the element grid associated with any available element.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, char *CODE, identifies the element desired.
    Legal values include "Q4", "Q8", "Q9", "Q12", "Q16", "QL", "T3", 
    "T4", "T6" and "T10".

    Input, int ELEMENT_ORDER, the number of nodes per element.

    Input, int NELEMX, NELEMY, the number of quadrilaterals along the
    X and Y directions.  The number of elements generated will be
    NELEMX * NELEMY for quadrilaterals, or 2 * NELEMX * NELEMY for
    triangles.

    Output, int GRID_ELEMENT[ELEMENT_ORDER*ELEMENT_NUM], the nodes 
    that form each element.  
*/
{
  int *element_node;

  if ( s_eqi ( code, "Q4" ) )
  {
    element_node = grid_q4_element ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "Q8" ) )
  {
    element_node = grid_q8_element ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "Q9" ) )
  {
    element_node = grid_q9_element ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "Q12" ) )
  {
    element_node = grid_q12_element ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "Q16" ) )
  {
    element_node = grid_q16_element ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "QL" ) )
  {
    element_node = grid_ql_element ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "T3" ) )
  {
    element_node = grid_t3_element ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "T4" ) )
  {
    element_node = grid_t4_element ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "T6" ) )
  {
    element_node = grid_t6_element ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "T10" ) )
  {
    element_node = grid_t10_element ( nelemx, nelemy );
  }
  else
  {
    element_node = NULL;
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "GRID_ELEMENT - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of CODE = \"%s\".\n", code );
    exit ( 1 );
  }

  return element_node;
}
/******************************************************************************/

int grid_element_num ( char *code, int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_ELEMENT_NUM returns the number of elements in a grid.

  Discussion:

    The number of elements generated will be NELEMX * NELEMY for
    quadrilaterals, or 2 * NELEMX * NELEMY for triangles.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, char *CODE, identifies the element desired.
    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 'T3', 
    'T4', 'T6' and 'T10'.

    Input, int NELEMX, NELEMY, the number of quadrilaterals along the
    X and Y directions.  

    Output, int GRID_ELEMENT_NUM, the number of elements in the grid.
*/
{
  int element_num;

  if ( s_eqi ( code, "Q4" ) )
  {
    element_num = grid_q4_element_num ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "Q8" ) )
  {
    element_num = grid_q8_element_num ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "Q9" ) )
  {
    element_num = grid_q9_element_num ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "Q12" ) )
  {
    element_num = grid_q12_element_num ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "Q16" ) )
  {
    element_num = grid_q16_element_num ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "QL" ) )
  {
    element_num = grid_ql_element_num ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "T3" ) )
  {
    element_num = grid_t3_element_num ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "T4" ) )
  {
    element_num = grid_t4_element_num ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "T6" ) )
  {
    element_num = grid_t6_element_num ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "T10" ) )
  {
    element_num = grid_t10_element_num ( nelemx, nelemy );
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "GRID_ELEMENT_NUM - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of CODE = \"%s\".\n", code );
    element_num = -1;
    exit ( 1 );
  }

  return element_num;
}
/******************************************************************************/

int grid_node_num ( char *code, int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_NODE_NUM returns the number of nodes in a grid.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, char *CODE, identifies the element desired.
    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 'T3', 
    'T4', 'T6' and 'T10'.

    Input, int NELEMX, NELEMY, the number of quadrilaterals along the
    X and Y directions.  

    Output, int GRID_NODE_NUM, the number of elements in the grid.
*/
{
  int node_num;

  if ( s_eqi ( code, "Q4" ) )
  {
    node_num = grid_q4_node_num ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "Q8" ) )
  {
    node_num = grid_q8_node_num ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "Q9" ) )
  {
    node_num = grid_q9_node_num ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "Q12" ) )
  {
    node_num = grid_q12_node_num ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "Q16" ) )
  {
    node_num = grid_q16_node_num ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "QL" ) )
  {
    node_num = grid_ql_node_num ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "T3" ) )
  {
    node_num = grid_t3_node_num ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "T4" ) )
  {
    node_num = grid_t4_node_num ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "T6" ) )
  {
    node_num = grid_t6_node_num ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "T10" ) )
  {
    node_num = grid_t10_node_num ( nelemx, nelemy );
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "GRID_NODE_NUM - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of CODE = \"%s\".\n", code );
    node_num = -1;
    exit ( 1 );
  }

  return node_num;
}
/******************************************************************************/

double *grid_nodes_01 ( int x_num, int y_num )

/******************************************************************************/
/*
  Purpose:

    GRID_NODES_01 returns an equally spaced rectangular grid of nodes in the unit square.

  Example:

    X_NUM = 5
    Y_NUM = 3

    NODE_XY =
    ( 0, 0.25, 0.5, 0.75, 1, 0,   0.25, 0.5, 0.75, 1,   0, 0.25, 0.5, 0.75, 1;
      0, 0,    0,   0,    0, 0.5, 0.5,  0.5, 0.5,  0.5, 1, 1.0,  1.0, 1.0,  1 )

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int X_NUM, Y_NUM, the number of nodes in the X and Y directions.

    Output, double GRID_NODES_01[2*X_NUM*Y_NUM], the coordinates of
    the nodes.
*/
{
  int i;
  int j;
  int k;
  int node_num;
  double *node_xy;
  double value;

  node_num = x_num * y_num;

  node_xy = ( double * ) malloc ( 2 * node_num * sizeof ( double ) );

  if ( x_num == 1 )
  {
    for ( k = 0; k < node_num; k++ )
    {
      node_xy[0+k*2] = 0.5;
    }
  }
  else
  {
    for ( i = 0; i < x_num; i++ )
    {
      value = ( double ) ( i ) / ( double ) ( x_num - 1 );
      for ( j = i; j < node_num; j = j + x_num )
      {
        node_xy[0+j*2] = value;
      }
    }
  }

  if ( y_num == 1 )
  {
    for ( k = 0; k < node_num; k++ )
    {
      node_xy[1+k*2] = 0.5;
    }
  }
  else
  {
    for ( j = 0; j < y_num; j++ )
    {
      value = ( double ) ( j ) / ( double ) ( y_num - 1 );
      for ( i = j*x_num; i < ( j + 1 ) * x_num; i++ )
      {
        node_xy[1+i*2] = value;
      }
    }
  }

  return node_xy;
}
/******************************************************************************/

void grid_print ( int element_order, int element_num, int element_node[] )

/******************************************************************************/
/*
  Purpose: 

    GRID_PRINT prints the elements that form a grid.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int ELEMENT_ORDER, the number of nodes per element.

    Input, int ELEMENT_NUM, the number of elements.
 
    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], the nodes that form
    each element.
*/
{
  int element;
  int i;

  printf ( "\n" );
  printf ( "  GRID_PRINT: Element -> Node table.\n" );
  printf ( "\n" );
  printf ( "  Element order = %d\n", element_order );
  printf ( "  Number of elements = %d\n", element_num );
  printf ( "\n" );
  printf ( "    #   " );
  for ( i = 0; i < element_order; i++ )
  {
    printf ( "%3d", i );
  }
  printf ( "\n" );
  printf ( "\n" );

  for ( element = 0; element < element_num; element++ )
  {
    printf ( "  %3d   ", element );
    for ( i = 0; i < element_order; i++ )
    {
      printf ( "%3d", element_node[i+element*element_order] );
    }
    printf ( "\n" );
  }

  return;
}
/******************************************************************************/

int *grid_q4_element ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose: 

    GRID_Q4_ELEMENT produces a grid of 4 node quadrilaterals.

  Discussion:

    For each element, the nodes are listed in counter-clockwise order.

  Example:

    Input:

      NELEMX = 3, NELEMY = 2

    Output:

      ELEMENT_NODE = 
         1, 2,  6,  5;
         2, 3,  7,  6;
         3, 4,  8,  7;
         5, 6, 10,  9;
         6, 7, 11, 10;
         7, 8, 12, 11.

  Grid:

    9---10---11---12
    |    |    |    |
    |    |    |    |
    |  4 |  5 |  6 |
    |    |    |    |
    5----6----7----8
    |    |    |    |
    |    |    |    |
    |  1 |  2 |  3 |
    |    |    |    |
    1----2----3----4

  Element Q4:

    |
    1  4-----3
    |  |     |
    |  |     |
    S  |     |
    |  |     |
    |  |     |
    0  1-----2
    |
    +--0--R--1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  The number of elements generated will be
    NELEMX * NELEMY.

    Output, int GRID_Q4[4*NELEMX*NELEMY], the nodes that form
    each element.
*/
{
  int element;
  int *element_node;
  int element_order = 4;
  int i;
  int j;
  int ne;
  int nw;
  int se;
  int sw;

  element_node = ( int * ) malloc ( element_order*nelemx*nelemy * sizeof ( int ) );
/*
  Node labeling:

    NW---NE
     |    |
    SW---SE
*/
  element = 0;
 
  for ( j = 1; j <= nelemy; j++ )
  {
    for ( i = 1; i <= nelemx; i++ )
    {
      sw = i     + ( j - 1 ) * ( nelemx + 1 );
      se = i + 1 + ( j - 1 ) * ( nelemx + 1 );
      nw = i     +   j       * ( nelemx + 1 );
      ne = i + 1 +   j       * ( nelemx + 1 );
  
      element_node[0+element*element_order] = sw;
      element_node[1+element*element_order] = se;
      element_node[2+element*element_order] = ne;
      element_node[3+element*element_order] = nw;

      element = element + 1;
    }
  }

  return element_node;
}
/******************************************************************************/

int grid_q4_element_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_Q4_ELEMENT_NUM counts the elements in a grid of 4 node quadrilaterals.

  Example:

    Input:

      NELEMX = 3, NELEMY = 2

    Output:

      ELEMENT_NUM = NELEMX * NELEMY = 6

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  The number of elements generated will be
    NELEMX * NELEMY.

    Output, int GRID_Q4_ELEMENT_NUM, the number of elements in the grid.
*/
{
  int element_num;

  element_num = nelemx * nelemy;

  return element_num;
}
/******************************************************************************/

int grid_q4_node_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_Q4_NODE_NUM counts the nodes in a grid of 4 node quadrilaterals.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  

    Output, int GRID_Q4_NODE_NUM, the number of nodes in the grid.
*/
{
  int node_num;

  node_num = ( nelemx + 1 ) * ( nelemy + 1 );

  return node_num;
}
/******************************************************************************/

int *grid_q8_element ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose: 

    GRID_Q8_ELEMENT produces a grid of 8 node quadrilaterals.

  Example:

    Input:

      NELEMX = 3, NELEMY = 2

    Output:

      ELEMENT_NODE =
         1,  3, 14, 12,  2,  9, 13,  8;
         3,  5, 16, 14,  4, 10, 15,  9;
         5,  7, 18, 16,  6, 11, 17, 10;
        12, 14, 25, 23, 13, 20, 24, 19;
        14, 16, 27, 25, 15, 21, 26, 20;
        16, 18, 29, 27, 17, 22, 28, 21.

  Diagram:

   23---24---25---26---27---28---29
    |         |         |         |
    |         |         |         |
   19        20        21        22
    |         |         |         |
    | 4       | 5       | 6       |
   12---13---14---15---16---17---18
    |         |         |         |
    |         |         |         |
    8         9        10        11
    |         |         |         |
    | 1       | 2       | 3       |
    1----2----3----4----5----6----7

  Element Q8:

    |
    1  4--7--3
    |  |     |
    |  |     |
    S  8     6
    |  |     |
    |  |     |
    0  1--5--2
    |
    +--0--R--1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  The number of elements generated will be
    NELEMX * NELEMY.

    Output, int GRID_Q8[8*NELEMX*NELEMY], the nodes that form
    each element.
*/
{
  int e;
  int element;
  int *element_node;
  int element_order = 8;
  int i;
  int j;
  int n;
  int ne;
  int nw;
  int s;
  int se;
  int sw;
  int w;

  element_node = ( int * ) malloc ( element_order*nelemx*nelemy * sizeof ( int ) );
/*
  Node labeling:

    NW----N----NE
     |          |
     W   (C)    E
     |          |
    SW----S----SE
*/

  element = 0;

  for ( j = 1; j <= nelemy; j++ )
  {
    for ( i = 1; i <= nelemx; i++ )
    {
      sw = ( j - 1 )  * ( 3 * nelemx + 2 ) + 2 * i - 1;
      w  = sw + 2 * nelemx + 2 - i;
      nw = sw + 3 * nelemx + 2;

      s =  sw + 1;
      n =  sw + ( 3 * nelemx + 2 ) + 1;

      se = sw + 2;
      e  = sw + 2 * nelemx + 2 - i + 1;
      ne = sw + ( 3 * nelemx + 2 ) + 2;

      element_node[0+element*element_order] = sw;
      element_node[1+element*element_order] = se;
      element_node[2+element*element_order] = ne;
      element_node[3+element*element_order] = nw;
      element_node[4+element*element_order] = s;
      element_node[5+element*element_order] = e;
      element_node[6+element*element_order] = n;
      element_node[7+element*element_order] = w;

      element = element + 1;
    }
  }

  return element_node;
}
/******************************************************************************/

int grid_q8_element_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_Q8_ELEMENT_NUM counts the elements in a grid of 8 node quadrilaterals.

  Example:

    Input:

      NELEMX = 3, NELEMY = 2

    Output:

      ELEMENT_NUM = NELEMX * NELEMY = 6

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  The number of elements generated will be
    NELEMX * NELEMY.

    Output, int GRID_Q8_ELEMENT_NUM, the number of elements in the grid.
*/
{
  int element_num;

  element_num = nelemx * nelemy;

  return element_num;
}
/******************************************************************************/

int grid_q8_node_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_Q8_NODE_NUM counts the nodes in a grid of 8 node quadrilaterals.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  

    Output, int GRID_Q8_NODE_NUM, the number of nodes in the grid.
*/
{
  int node_num;

  node_num = 3 * nelemx * nelemy + 2 * nelemx + 2 * nelemy + 1;

  return node_num;
}
/******************************************************************************/

int *grid_q9_element ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose: 

    GRID_Q9_ELEMENT produces a grid of 9 node quadrilaterals.

  Example:

    Input:

      NELEMX = 3, NELEMY = 2

    Output:

      ELEMENT_NODE = 
         1,  3, 17, 15,  2, 10, 16,  8,  9;
         3,  5, 19, 17,  4, 12, 18, 10, 11;
         5,  7, 21, 19,  6, 14, 20, 12, 13;
        15, 17, 31, 29, 16, 24, 30, 22, 23;
        17, 19, 33, 31, 18, 26, 32, 24, 25;
        19, 21, 35, 33, 20, 28, 34, 26, 27.

  Grid:

   29---30---31---32---33---34---35
    |    .    |    .    |    .    |
    |    .    |    .    |    .    |
   22 . 23 . 24 . 25 . 26 . 27 . 28
    |    .    |    .    |    .    |
    | 4  .    | 5  .    | 6  .    |
   15---16---17---18---19---20---21
    |    .    |    .    |    .    |
    |    .    |    .    |    .    |
    8 .  9 . 10 . 11 . 12 . 13 . 14
    |    .    |    .    |    .    |
    | 1  .    | 2  .    | 3  .    |
    1----2----3----4----5----6----7

  Element Q9:

    |
    1  4--7--3
    |  |     |
    |  |     |
    S  8  9  6
    |  |     |
    |  |     |
    0  1--5--2
    |
    +--0--R--1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  The number of elements generated will be
    NELEMX * NELEMY.

    Output, int GRID_Q9[9*NELEMX*NELEMY], the nodes that form
    each element.
*/
{
  int c;
  int e;
  int element;
  int *element_node;
  int element_order = 9;
  int i;
  int j;
  int n;
  int ne;
  int nw;
  int s;
  int se;
  int sw;
  int w;

  element_node = ( int * ) malloc ( element_order*nelemx*nelemy * sizeof ( int ) );
/*
  Node labeling:

    NW----N----NE
     |          |
     W    C     E
     |          |
    SW----S----SE
*/
  element = 0;

  for ( j = 1; j <= nelemy; j++ )
  {
    for ( i = 1; i <= nelemx; i++ )
    {
      sw = 2 * ( j - 1 )  * ( 2 * nelemx + 1 ) + 2 * ( i - 1 ) + 1;
      w  = sw +               2 * nelemx + 1;
      nw = sw +         2 * ( 2 * nelemx + 1 );

      s  = sw + 1;
      c  = sw + 1 +               2 * nelemx + 1;
      n  = sw + 1 +         2 * ( 2 * nelemx + 1 );

      se = sw + 2;
      e  = sw + 2 +               2 * nelemx + 1;
      ne = sw + 2 +         2 * ( 2 * nelemx + 1 );

      element_node[0+element*element_order] = sw;
      element_node[1+element*element_order] = se;
      element_node[2+element*element_order] = ne;
      element_node[3+element*element_order] = nw;
      element_node[4+element*element_order] = s;
      element_node[5+element*element_order] = e;
      element_node[6+element*element_order] = n;
      element_node[7+element*element_order] = w;
      element_node[8+element*element_order] = c;

      element = element + 1;
    }
  }

  return element_node;
}
/******************************************************************************/

int grid_q9_element_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_Q9_ELEMENT_NUM counts the elements in a grid of 9 node quadrilaterals.

  Example:

    Input:

      NELEMX = 3, NELEMY = 2

    Output:

      ELEMENT_NUM = NELEMX * NELEMY = 6

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  The number of elements generated will be
    NELEMX * NELEMY.

    Output, int GRID_Q9_ELEMENT_NUM, the number of elements in the grid.
*/
{
  int element_num;

  element_num = nelemx * nelemy;

  return element_num;
}
/******************************************************************************/

int grid_q9_node_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_Q9_NODE_NUM counts the nodes in a grid of 9 node quadrilaterals.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    23 March 2006

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  

    Output, int GRID_Q9_NODE_NUM, the number of nodes in the grid.
*/
{
  int node_num;

  node_num = ( 2 * nelemx + 1 ) * ( 2 * nelemy + 1 );

  return node_num;
}
/******************************************************************************/

int *grid_q12_element ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose: 

    GRID_Q12_ELEMENT produces a grid of 12 node quadrilaterals.

  Example:

    Input:

      NELEMX = 3, NELEMY = 2

    Output:

      ELEMENT_NODE =
         1,  2,  3,  4, 11, 12, 15, 16, 19, 20, 21, 22;
         4,  5,  6,  7, 12, 13, 16, 17, 22, 23, 24, 25;
         7,  8,  9, 10, 13, 14, 17, 18, 25, 26, 27, 28;
        19, 20, 21, 22, 29, 30, 33, 34, 37, 38, 39, 40;
        22, 23, 24, 25, 30, 31, 34, 35, 40, 41, 42, 43;
        25, 26, 27, 28, 31, 32, 35, 36, 43, 44, 45, 46.

  Grid:

   37-38-39-40-41-42-43-44-45-46
    |        |        |        |
   33       34       35       36
    |        |        |        |
   29       30       31       32
    | 4      | 5      | 6      |
   19-20-21-22-23-24-25-26-27-28
    |        |        |        |
   15       16       17       18
    |        |        |        |
   11       12       13       14
    | 1      | 2      | 3      |
    1--2--3--4--5--6--7--8--9-10

  Element Q12:

    |
    1  9-10-11-12
    |  |        |
    |  7        8
    S  |        |
    |  5        6
    |  |        |
    0  1--2--3--4
    |
    +--0---R---1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  The number of elements generated will be
    NELEMX * NELEMY.

    Output, int GRID_Q12[12*NELEMX*NELEMY], the nodes that form
    each element.
*/
{
  int base;
  int element;
  int *element_node;
  int element_order = 12;
  int i;
  int j;

  element_node = ( int * ) malloc ( element_order*nelemx*nelemy * sizeof ( int ) );

  element = 0;

  for ( j = 1; j <= nelemy; j++ )
  {
    for ( i = 1; i <= nelemx; i++ )
    {
      base = ( j - 1 )  * ( 5 * nelemx + 3 ) + 1;

      element_node[ 0+element*element_order] =  base + ( i - 1 ) * 3;
      element_node[ 1+element*element_order] =  base + ( i - 1 ) * 3 + 1;
      element_node[ 2+element*element_order] =  base + ( i - 1 ) * 3 + 2;
      element_node[ 3+element*element_order] =  base + ( i - 1 ) * 3 + 3;

      element_node[ 4+element*element_order] =  base + 3 * nelemx + i;
      element_node[ 5+element*element_order] =  base + 3 * nelemx + i + 1;

      element_node[ 6+element*element_order] =  base + 4 * nelemx + i + 1;
      element_node[ 7+element*element_order] =  base + 4 * nelemx + i + 2;

      element_node[ 8+element*element_order] =  base + 5 * nelemx + 3 * i;
      element_node[ 9+element*element_order] = base + 5 * nelemx + 3 * i + 1;
      element_node[10+element*element_order] = base + 5 * nelemx + 3 * i + 2;
      element_node[11+element*element_order] = base + 5 * nelemx + 3 * i + 3;

      element = element + 1;
    }
  }

  return element_node;
}
/******************************************************************************/

int grid_q12_element_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_Q12_ELEMENT_NUM counts the elements in a grid of 12 node quadrilaterals.

  Example:

    Input:

      NELEMX = 3, NELEMY = 2

    Output:

      ELEMENT_NUM = NELEMX * NELEMY = 6

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  The number of elements generated will be
    NELEMX * NELEMY.

    Output, int GRID_Q12_ELEMENT_NUM, the number of elements in the grid.
*/
{
  int element_num;

  element_num = nelemx * nelemy;

  return element_num;
}
/******************************************************************************/

int grid_q12_node_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_Q12_NODE_NUM counts the nodes in a grid of 12 node quadrilaterals.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  

    Output, int GRID_Q12_NODE_NUM, the number of nodes in the grid.
*/
{
  int node_num;

  node_num = 5 * nelemx * nelemy + 3 * nelemx + 3 * nelemy + 1;

  return node_num;
}
/******************************************************************************/

int *grid_q16_element ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose: 

    GRID_Q16_ELEMENT produces a grid of 16 node quadrilaterals.

  Example:

    Input:

      NELEMX = 2, NELEMY = 2

    Output:

      ELEMENT_NODE = 
         1,  2,  3,  4,  8,  9, 10, 11, 15, 16, 17, 18, 22, 23, 24, 25;
         4,  5,  6,  7, 11, 12, 13, 14, 18, 19, 20, 21, 25, 26, 27, 28;
        22, 23, 24, 25, 29, 30, 31, 32, 36, 37, 38, 39, 43, 44, 45, 46;
        25, 26, 27, 28, 32, 33, 34, 35, 39, 40, 41, 42, 46, 47, 48, 49. 
        
  Grid:

   43-44-45-46-47-48-49
    |        |        |
    |        |        |
   36 37 38 39 40 41 42
    |        |        |
    |        |        |
   29 30 31 32 33 34 35
    |        |        |
    | 3      | 4      |
   22-23-24-25-26-27-28
    |        |        |
    |        |        |
   15 16 17 18 19 20 21
    |        |        |
    |        |        |
    8  9 10 11 12 13 14
    |        |        |
    | 1      | 2      |
    1--2--3--4--5--6--7

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  The number of elements generated will be
    NELEMX * NELEMY.

    Output, int GRID_Q16[16*NELEMX*NELEMY], the nodes that form
    each element.
*/
{
  int base;
  int element;
  int *element_node;
  int element_order = 16;
  int i;
  int j;

  element_node = ( int * ) malloc ( element_order*nelemx*nelemy * sizeof ( int ) );

  element = 0;
 
  for ( j = 1; j <= nelemy; j++ )
  {
    for ( i = 1; i <= nelemx; i++ )
    {
      base = ( j - 1 ) * 3 * ( 3 * nelemx + 1 ) + 3 * i - 2;

      element_node[ 0+element*element_order] = base;
      element_node[ 1+element*element_order] = base                          + 1;
      element_node[ 2+element*element_order] = base                          + 2;
      element_node[ 3+element*element_order] = base                          + 3;
      element_node[ 4+element*element_order] = base +     ( 3 * nelemx + 1 );
      element_node[ 5+element*element_order] = base +     ( 3 * nelemx + 1 ) + 1;
      element_node[ 6+element*element_order] = base +     ( 3 * nelemx + 1 ) + 2;
      element_node[ 7+element*element_order] = base +     ( 3 * nelemx + 1 ) + 3;
      element_node[ 8+element*element_order] = base + 2 * ( 3 * nelemx + 1 );
      element_node[ 9+element*element_order] = base + 2 * ( 3 * nelemx + 1 ) + 1;
      element_node[10+element*element_order] = base + 2 * ( 3 * nelemx + 1 ) + 2;
      element_node[11+element*element_order] = base + 2 * ( 3 * nelemx + 1 ) + 3;
      element_node[12+element*element_order] = base + 3 * ( 3 * nelemx + 1 );
      element_node[13+element*element_order] = base + 3 * ( 3 * nelemx + 1 ) + 1;
      element_node[14+element*element_order] = base + 3 * ( 3 * nelemx + 1 ) + 2;
      element_node[15+element*element_order] = base + 3 * ( 3 * nelemx + 1 ) + 3;

      element = element + 1;
    }
  }

  return element_node;
}
/******************************************************************************/

int grid_q16_element_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_Q16_ELEMENT_NUM counts the elements in a grid of 16 node quadrilaterals.

  Example:

    Input:

      NELEMX = 3, NELEMY = 2

    Output:

      ELEMENT_NUM = NELEMX * NELEMY = 6

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  The number of elements generated will be
    NELEMX * NELEMY.

    Output, int GRID_Q16_ELEMENT_NUM, the number of elements in the grid.
*/
{
  int element_num;

  element_num = nelemx * nelemy;

  return element_num;
}
/******************************************************************************/

int grid_q16_node_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_Q16_NODE_NUM counts the nodes in a grid of 16 node quadrilaterals.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  

    Output, int GRID_Q16_NODE_NUM, the number of nodes in the grid.
*/
{
  int node_num;

  node_num = ( 3 * nelemx + 1 ) * ( 3 * nelemy + 1 );

  return node_num;
}
/******************************************************************************/

int *grid_ql_element ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose: 

    GRID_QL_ELEMENT produces a grid of 6 node quadratics/linears.

  Example:

    Input:

      NELEMX = 3, NELEMY = 2

    Output:

      ELEMENT_NODE = 
         1,  2,  3,  8,  9, 10;
         3,  4,  5, 10, 11, 12;
         5,  6,  7, 12, 13, 14;
         8,  9, 10, 15, 16, 17;
        10, 11, 12, 17, 18, 19;
        12, 13, 14, 19, 20, 21.

  Grid:

   15---16---17---18---19---20---21
    |         |         |         |
    |         |         |         |
    |    4    |    5    |    6    |
    |         |         |         |
    |         |         |         |
    8----9---10---11---12---13---14
    |         |         |         |
    |         |         |         |
    |    1    |    2    |    3    |
    |         |         |         |
    |         |         |         |
    1----2----3----4----5----6----7

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  The number of elements generated will be
    NELEMX * NELEMY.  X will the the "quadratic direction", and
    Y will be the "linear direction".

    Output, int GRID_QL[6*NELEMX*NELEMY], the nodes that form
    each element.
*/
{
  int base;
  int element;
  int *element_node;
  int element_order = 6;
  int i;
  int j;

  element_node = ( int * ) malloc ( element_order*nelemx*nelemy * sizeof ( int ) );

  element = 0;
 
  for ( j = 1; j <= nelemy; j++ )
  {
    for ( i = 1; i <= nelemx; i++ )
    {
      base = ( j - 1 )  * ( 2 * nelemx + 1 ) + 2 * i - 1;

      element_node[0+element*element_order] = base;
      element_node[1+element*element_order] = base + 1;
      element_node[2+element*element_order] = base + 2;
      element_node[3+element*element_order] = base + ( 2 * nelemx + 1 );
      element_node[4+element*element_order] = base + ( 2 * nelemx + 1 ) + 1;
      element_node[5+element*element_order] = base + ( 2 * nelemx + 1 ) + 2;

      element = element + 1;
    }
  }

  return element_node;
}
/******************************************************************************/

int grid_ql_element_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_QL_ELEMENT_NUM counts the elements in a grid of quadratic/linear quadrilaterals.

  Example:

    Input:

      NELEMX = 3, NELEMY = 2

    Output:

      ELEMENT_NUM = NELEMX * NELEMY = 6

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  The number of elements generated will be
    NELEMX * NELEMY.

    Output, int GRID_QL_ELEMENT_NUM, the number of elements in the grid.
*/
{
  int element_num;

  element_num = nelemx * nelemy;

  return element_num;
}
/******************************************************************************/

int grid_ql_node_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_QL_NODE_NUM counts the nodes in a grid of quadratic/linear quadrilaterals.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  

    Output, int GRID_QL_NODE_NUM, the number of nodes in the grid.
*/
{
  int node_num;

  node_num = 2 * nelemx * nelemy + 2 * nelemx + nelemy + 1;

  return node_num;
}
/******************************************************************************/

void grid_shape_2d ( int n, double a[], int *n1, int *n2 )

/******************************************************************************/
/*
  Purpose: 

    GRID_SHAPE_2D guesses the shape N1 by N2 of a vector of data.

  Discussion:

    The data vector A is assumed to contain N1 * N2 values, with
    where each of N2 values is repeated N1 times.

  Example:

    Input:

      A = ( 2, 2, 2, 7, 7, 7 )

    Output:

      N1 = 3, N2 = 2

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of data values.

    Input, double A[N], the data, which should have the properties
    described above.

    Output, int *N1, *N2, the "shape" of the data in the array.
*/
{
  int i;
/*
  Make a guess for N1.
*/
  i = 1;
  *n1 = 1;

  for ( i = 1; i < n; i++ )
  {
    if ( a[i] != a[0] )
    {
      break;
    }
    *n1 = *n1 + 1;
  }
/*
  Guess that N2 = N / N1.
*/
  *n2 = n / (*n1);

  return;
}
/******************************************************************************/

int *grid_t3_element ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose: 

    GRID_T3_ELEMENT produces a grid of pairs of 3 node triangles.

  Example:

    Input:

      NELEMX = 3, NELEMY = 2

    Output:

      ELEMENT_NODE = 
         1,  2,  5;
         6,  5,  2;
         2,  3,  6;
         7,  6,  3;
         3,  4,  7;
         8,  7,  4;
         5,  6,  9;
        10,  9,  6;
         6,  7, 10;
        11, 10,  7;
         7,  8, 11;
        12, 11,  8.

  Grid:

    9---10---11---12
    |\ 8 |\10 |\12 |
    | \  | \  | \  |
    |  \ |  \ |  \ |
    |  7\|  9\| 11\|
    5----6----7----8
    |\ 2 |\ 4 |\ 6 |
    | \  | \  | \  |
    |  \ |  \ |  \ |
    |  1\|  3\|  5\|
    1----2----3----4

  Element T3:

    |
    1  3
    |  |\
    |  | \
    S  |  \
    |  |   \
    |  |    \
    0  1-----2
    |
    +--0--R--1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  The number of elements generated will be
    2 * NELEMX * NELEMY.

    Output, int GRID_T3[3*2*NELEMX*NELEMY], the nodes that form
    each element.
*/
{
  int element;
  int *element_node;
  int element_order = 3;
  int i;
  int j;
  int ne;
  int nw;
  int se;
  int sw;

  element_node = ( int * ) malloc ( element_order*2*nelemx*nelemy * sizeof ( int ) );
/*
  Node labeling:

    NW--NE
     |\ |
     | \|
    SW--SE
*/
  element = 0;
 
  for ( j = 1; j <= nelemy; j++ )
  {
    for ( i = 1; i <= nelemx; i++ )
    {
      sw = i     + ( j - 1 ) * ( nelemx + 1 );
      se = i + 1 + ( j - 1 ) * ( nelemx + 1 );
      nw = i     +   j       * ( nelemx + 1 );
      ne = i + 1 +   j       * ( nelemx + 1 );

      element_node[0+element*element_order] = sw;
      element_node[1+element*element_order] = se;
      element_node[2+element*element_order] = nw;
      element = element + 1;

      element_node[0+element*element_order] = ne;
      element_node[1+element*element_order] = nw;
      element_node[2+element*element_order] = se;
      element = element + 1;
    }
  }

  return element_node;
}
/******************************************************************************/

int grid_t3_element_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_T3_ELEMENT_NUM counts the elements in a grid of 3 node triangles.

  Example:

    Input:

      NELEMX = 3, NELEMY = 2

    Output:

      ELEMENT_NUM = 2 * NELEMX * NELEMY = 12

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  The number of elements generated will be
    2 * NELEMX * NELEMY.

    Output, int GRID_T3_ELEMENT_NUM, the number of elements in the grid.
*/
{
  int element_num;

  element_num = 2 * nelemx * nelemy;

  return element_num;
}
/******************************************************************************/

int grid_t3_node_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_T3_NODE_NUM counts the nodes in a grid of 3 node triangles.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  

    Output, int GRID_T3_NODE_NUM, the number of nodes in the grid.
*/
{
  int node_num;

  node_num = ( nelemx + 1 ) * ( nelemy + 1 );

  return node_num;
}
/******************************************************************************/

int *grid_t4_element ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose: 

    GRID_T4_ELEMENT produces a grid of pairs of 4 node triangles.

  Example:

    Input:

      NELEMX = 3, NELEMY = 2

    Output:

      ELEMENT_NODE = 
         1,  2,  11,  5;
        12, 11,   2,  8;
         2,  3,  12,  6;
        13, 12,   3,  9;
         3   4   13,  7;
        14, 13,   4,  10;
        11, 12,  21,  15;
        22, 21,  12,  18;
        12, 13,  22,  16;
        23, 22,  13,  19;
        13  14   23,  17;
        24, 23,  14,  20;

  Grid:

   21---22---23---24
    |\18 |\19 |\20 |
    | \  | \  | \  |
    |  \ |  \ |  \ |
    | 15\| 16\| 17\|
   11---12---13---14
    |\ 8 |\ 9 |\10 |
    | \  | \  | \  |
    |  \ |  \ |  \ |
    | 5 \|  6\|  7\|
    1----2----3----4

  Element T4:

    |
    1  3
    |  |\
    |  | \
    S  |  \
    |  | 4 \
    |  |    \
    0  1-----2
    |
    +--0--R--1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  The number of elements generated will be
    2 * NELEMX * NELEMY.

    Output, int GRID_T4[4*2*NELEMX*NELEMY], the nodes that form
    each element.
*/
{
  int element;
  int *element_node;
  int element_order = 4;
  int i;
  int j;
  int nc;
  int ne;
  int nw;
  int sc;
  int se;
  int sw;

  element_node = ( int * ) malloc ( element_order*2*nelemx*nelemy * sizeof ( int ) );
/*
  Node labeling:

    NW----NE
     |\   |
     | \NC|
     |SC\ |
     |   \|
    SW---SE
*/
  element = 0;
 
  for ( j = 1; j <= nelemy; j++ )
  {
    for ( i = 1; i <= nelemx; i++ )
    {
      sw = i     + ( j - 1 ) * ( 3 * nelemx + 1 );
      se = sw + 1;
      sc = sw +     nelemx + 1;
      nc = sw + 2 * nelemx + 1;
      nw = sw + 3 * nelemx + 1;
      ne = sw + 3 * nelemx + 2;

      element_node[0+element*element_order] = sw;
      element_node[1+element*element_order] = se;
      element_node[2+element*element_order] = nw;
      element_node[3+element*element_order] = sc;
      element = element + 1;

      element_node[0+element*element_order] = ne;
      element_node[1+element*element_order] = nw;
      element_node[2+element*element_order] = se;
      element_node[3+element*element_order] = nc;
      element = element + 1;
    }
  }

  return element_node;
}
/******************************************************************************/

int grid_t4_element_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_T4_ELEMENT_NUM counts the elements in a grid of 4 node triangles.

  Example:

    Input:

      NELEMX = 3, NELEMY = 2

    Output:

      ELEMENT_NUM = 2 * NELEMX * NELEMY = 12

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  The number of elements generated will be
    2 * NELEMX * NELEMY.

    Output, int GRID_T4_ELEMENT_NUM, the number of elements in the grid.
*/
{
  int element_num;

  element_num = 2 * nelemx * nelemy;

  return element_num;
}
/******************************************************************************/

int grid_t4_node_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_T4_NODE_NUM counts the nodes in a grid of 4 node triangles.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  

    Output, int GRID_T4_NODE_NUM, the number of nodes in the grid.
*/
{
  int node_num;

  node_num = ( nelemx + 1 ) * ( nelemy + 1 ) + 2 * nelemx * nelemy;

  return node_num;
}
/******************************************************************************/

int *grid_t6_element ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose: 

    GRID_T6_ELEMENT produces a grid of pairs of 6 node triangles.

  Example:

    Input:

      NELEMX = 3, NELEMY = 2

    Output:

      ELEMENT_NODE = 
         1,  3, 15,  2,  9,  8;
        17, 15,  3, 16,  9, 10;
         3,  5, 17,  4, 11, 10;
        19, 17,  5, 18, 11, 12;
         5,  7, 19,  6, 13, 12;
        21, 19,  7, 20, 13, 14;
        15, 17, 29, 16, 23, 22;
        31, 29, 17, 30, 23, 24;
        17, 19, 31, 18, 25, 24;
        33, 31, 19, 32, 25, 26;
        19, 21, 33, 20, 27, 26;
        35, 33, 21, 34, 27, 28.

  Grid:

   29-30-31-32-33-34-35
    |\ 8  |\10  |\12  |
    | \   | \   | \   |
   22 23 24 25 26 27 28
    |   \ |   \ |   \ |
    |  7 \|  9 \| 11 \|
   15-16-17-18-19-20-21
    |\ 2  |\ 4  |\ 6  |
    | \   | \   | \   |
    8  9 10 11 12 13 14
    |   \ |   \ |   \ |
    |  1 \|  3 \|  5 \|
    1--2--3--4--5--6--7

  Element T6:

    |
    1  3
    |  |\
    |  | \
    S  6  5
    |  |   \
    |  |    \
    0  1--4--2
    |
    +--0--R--1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  The number of elements generated will be
    2 * NELEMX * NELEMY.

    Output, int GRID_T6[6*2*NELEMX*NELEMY], the nodes that form
    each element.
*/
{
  int c;
  int e;
  int element;
  int *element_node;
  int element_order = 6;
  int i;
  int j;
  int n;
  int ne;
  int nw;
  int s;
  int se;
  int sw;
  int w;

  element_node = ( int * ) malloc ( element_order*2*nelemx*nelemy * sizeof ( int ) );
/*
  Node labeling:

    NW---N--NE
     | \     |
     W   C   E
     |    \  |
    SW---S--SE
*/
  element = 0;
 
  for ( j = 1; j <= nelemy; j++ )
  {
    for ( i = 1; i <= nelemx; i++ )
    {
      sw = 2 * ( j - 1 )  * ( 2 * nelemx + 1 ) + 2 * ( i - 1 ) + 1;
      w  = sw +               2 * nelemx + 1;
      nw = sw +         2 * ( 2 * nelemx + 1 );

      s  = sw + 1;
      c  = sw + 1 +               2 * nelemx + 1;
      n  = sw + 1 +         2 * ( 2 * nelemx + 1 );

      se = sw + 2;
      e  = sw + 2 +               2 * nelemx + 1;
      ne = sw + 2 +         2 * ( 2 * nelemx + 1 );

      element_node[0+element*element_order] = sw;
      element_node[1+element*element_order] = se;
      element_node[2+element*element_order] = nw;
      element_node[3+element*element_order] = s;
      element_node[4+element*element_order] = c;
      element_node[5+element*element_order] = w;
      element = element + 1;

      element_node[0+element*element_order] = ne;
      element_node[1+element*element_order] = nw;
      element_node[2+element*element_order] = se;
      element_node[3+element*element_order] = n;
      element_node[4+element*element_order] = c;
      element_node[5+element*element_order] = e;
      element = element + 1;
    }
  }

  return element_node;
}
/******************************************************************************/

int grid_t6_element_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_T6_ELEMENT_NUM counts the elements in a grid of 6 node triangles.

  Example:

    Input:

      NELEMX = 3, NELEMY = 2

    Output:

      ELEMENT_NUM = 2 * NELEMX * NELEMY = 12

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  The number of elements generated will be
    2 * NELEMX * NELEMY.

    Output, int GRID_T6_ELEMENT_NUM, the number of elements in the grid.
*/
{
  int element_num;

  element_num = 2 * nelemx * nelemy;

  return element_num;
}
/******************************************************************************/

int grid_t6_node_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_T6_NODE_NUM counts the nodes in a grid of 6 node triangles.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  

    Output, int GRID_T6_NODE_NUM, the number of nodes in the grid.
*/
{
  int node_num;

  node_num = ( 2 * nelemx + 1 ) * ( 2 * nelemy + 1 );

  return node_num;
}
/******************************************************************************/

int *grid_t10_element ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose: 

    GRID_T10_ELEMENT produces a grid of pairs of 10 node triangles.

  Example:

    Input:

      NELEMX = 2, NELEMY = 2

    Output:

      ELEMENT_NODE = 
         1,  2,  3,  4, 10, 16, 22, 15,  8,  9;
        25, 24, 23, 22, 16, 10,  4, 11, 18, 17;
         4,  5,  6,  7, 13, 19, 25, 18, 11, 12;
        28, 27, 26, 25, 19, 13,  7, 14, 21, 20;
        22, 23, 24, 25, 31, 37, 43, 36, 29, 30;
        46, 45, 44, 43, 37, 31, 25, 32, 39, 38;
        25, 26, 27, 28, 34, 40, 46, 39, 31, 33;
        49, 48, 47, 46, 40, 34, 28, 35, 42, 41.
        
  Grid:

   43-44-45-46-47-48-49
    |\     6 |\     8 |
    | \      | \      |
   36 37 38 39 40 41 42
    |   \    |   \    |
    |    \   |    \   |
   29 30 31 32 33 34 35
    |      \ |      \ |
    | 5     \| 7     \|
   22-23-24-25-26-27-28
    |\     2 |\     4 |
    | \      | \      |
   15 16 17 18 19 20 21
    |   \    |   \    |
    |    \   |    \   |
    8  9 10 11 12 13 14
    |      \ |      \ |
    | 1     \| 3     \|
    1--2--3--4--5--6--7

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  The number of elements generated will be
    2 * NELEMX * NELEMY.

    Output, int GRID_T10[10*2*NELEMX*NELEMY], the nodes that form
    each element.
*/
{
  int base;
  int element;
  int *element_node;
  int element_order = 10;
  int i;
  int j;

  element_node = ( int * ) malloc ( element_order*2*nelemx*nelemy * sizeof ( int ) );

  element = 0;
 
  for ( j = 1; j <= nelemy; j++ )
  {
    for ( i = 1; i <= nelemx; i++ )
    {
      base = ( j - 1 ) * 3 * ( 3 * nelemx + 1 ) + 3 * i - 2;

      element_node[0+element*element_order] = base;
      element_node[1+element*element_order] = base                          + 1;
      element_node[2+element*element_order] = base                          + 2;
      element_node[3+element*element_order] = base                          + 3;
      element_node[4+element*element_order] = base +     ( 3 * nelemx + 1 ) + 2;
      element_node[5+element*element_order] = base + 2 * ( 3 * nelemx + 1 ) + 1;
      element_node[6+element*element_order] = base + 3 * ( 3 * nelemx + 1 );
      element_node[7+element*element_order] = base + 2 * ( 3 * nelemx + 1 );
      element_node[8+element*element_order] = base +     ( 2 * nelemx + 1 ) + 2;
      element_node[9+element*element_order] = base +     ( 2 * nelemx + 1 ) + 3;
      element = element + 1;

      element_node[0+element*element_order] = base + 3 * ( 3 * nelemx + 1 ) + 3;
      element_node[1+element*element_order] = base + 3 * ( 3 * nelemx + 1 ) + 2;
      element_node[2+element*element_order] = base + 3 * ( 3 * nelemx + 1 ) + 1;
      element_node[3+element*element_order] = base + 3 * ( 3 * nelemx + 1 );
      element_node[4+element*element_order] = base + 2 * ( 3 * nelemx + 1 ) + 1;
      element_node[5+element*element_order] = base +     ( 3 * nelemx + 1 ) + 2;
      element_node[6+element*element_order] = base                          + 3;
      element_node[7+element*element_order] = base +     ( 3 * nelemx + 1 ) + 3;
      element_node[8+element*element_order] = base + 2 * ( 3 * nelemx + 1 ) + 3;
      element_node[9+element*element_order] = base + 2 * ( 3 * nelemx + 1 ) + 2;
      element = element + 1;
    }
  }

  return element_node;
}
/******************************************************************************/

int grid_t10_element_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_T10_ELEMENT_NUM counts the elements in a grid of 10 node triangles.

  Example:

    Input:

      NELEMX = 3, NELEMY = 2

    Output:

      ELEMENT_NUM = 2 * NELEMX * NELEMY = 12

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  The number of elements generated will be
    2 * NELEMX * NELEMY.

    Output, int GRID_T10_ELEMENT_NUM, the number of elements in the grid.
*/
{
  int element_num;

  element_num = 2 * nelemx * nelemy;

  return element_num;
}
/******************************************************************************/

int grid_t10_node_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_T10_NODE_NUM counts the nodes in a grid of 10 node triangles.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  

    Output, int GRID_T10_NODE_NUM, the number of nodes in the grid.
*/
{
  int node_num;

  node_num = ( 3 * nelemx + 1 ) * ( 3 * nelemy + 1 );

  return node_num;
}
/******************************************************************************/

void grid_test ( char *code )

/******************************************************************************/
/*
  Purpose: 

    GRID_TEST tests the grid routines.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 July 2013

  Author:

    John Burkardt

  Parameters:

    Input, char *CODE, the code for the element.
    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 
    'T3', 'T4', 'T6' and 'T10'.
*/
{
  int *element_node;
  int element_num;
  int element_order;
  int nelemx;
  int nelemy;
  int width;
/*
  NODE is defined as a vector rather than a two dimensional array,
  so that we can handle the various cases using a single array.
*/
  printf ( "\n" );
  printf ( "  GRID_TEST: Test the grid routine for element %s\n", code );

  nelemx = 3;
  nelemy = 2;

  if ( s_eqi ( code, "Q4" ) || 
       s_eqi ( code, "Q8" ) || 
       s_eqi ( code, "Q9" ) ||
       s_eqi ( code, "Q12" ) ||
       s_eqi ( code, "Q16" ) ||
       s_eqi ( code, "QL" ) )
  {
    element_num = nelemx * nelemy;
  }
  else if ( strcmp ( code, "T3" ) == 0 || 
            strcmp ( code, "T4" ) == 0 || 
            strcmp ( code, "T6" ) == 0 || 
            strcmp ( code, "T10" ) == 0 )
  {
    element_num = 2 * nelemx * nelemy;
  }

  element_order = order_code ( code );

  element_node = grid_element ( code, element_order, nelemx, nelemy );

  grid_print ( element_order, element_num, element_node );

  width = grid_width ( element_order, element_num, element_node );

  printf ( "\n" );
  printf ( "  Grid width is %d\n", width );

  free ( element_node );

  return;
}
/******************************************************************************/

int grid_width ( int element_order, int element_num, int element_node[] )

/******************************************************************************/
/*
  Purpose: 

    GRID_WIDTH computes the width of a given grid.

  Definition:

    The grid width is defined to the maximum absolute
    difference of global indices of nodes in the same element.

  Example:

    For the following grid, the grid width is 13.

   23---24---25---26---27---28---29
    |         |         |         |
    |         |         |         |
   19        20        21        22
    |         |         |         |
    | 4       | 5       | 6       |
   12---13---14---15---16---17---18
    |         |         |         |
    |         |         |         |
    8         9        10        11
    |         |         |         |
    | 1       | 2       | 3       |
    1----2----3----4----5----6----7

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int ELEMENT_ORDER, the order of the elements.

    Input, int ELEMENT_NUM, the number of elements.
 
    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], the nodes that form
    each element.

    Output, int GRID_WIDTH, the grid width.
*/
{
  int element;
  int ip1;
  int ip2;
  int node1;
  int node2;
  int width;

  width = 0;
 
  for ( element = 0; element < element_num; element++ )
  {
    for ( node1 = 0; node1 < element_order; node1++ )
    {
      ip1 = element_node[node1+element*element_order];
      for ( node2 = 0; node2 < element_order; node2++ )
      {
        ip2 = element_node[node2+element*element_order];
        width = i4_max ( width, abs ( ip1 - ip2 ) );
      }
    }
  }
 
  return width;
}
/******************************************************************************/

int i4_max ( int i1, int i2 )

/******************************************************************************/
/*
  Purpose:

    I4_MAX returns the maximum of two I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 August 2006

  Author:

    John Burkardt

  Parameters:

    Input, int I1, I2, are two integers to be compared.

    Output, int I4_MAX, the larger of I1 and I2.
*/
{
  int value;

  if ( i2 < i1 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
/******************************************************************************/

int i4_min ( int i1, int i2 )

/******************************************************************************/
/*
  Purpose:

    I4_MIN returns the smaller of two I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 August 2006

  Author:

    John Burkardt

  Parameters:

    Input, int I1, I2, two integers to be compared.

    Output, int I4_MIN, the smaller of I1 and I2.
*/
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
/******************************************************************************/

int i4_modp ( int i, int j )

/******************************************************************************/
/*
  Purpose:

    I4_MODP returns the nonnegative remainder of I4 division.

  Discussion:

    If
      NREM = I4_MODP ( I, J )
      NMULT = ( I - NREM ) / J
    then
      I = J * NMULT + NREM
    where NREM is always nonnegative.

    The MOD function computes a result with the same sign as the
    quantity being divided.  Thus, suppose you had an angle A,
    and you wanted to ensure that it was between 0 and 360.
    Then mod(A,360) would do, if A was positive, but if A
    was negative, your result would be between -360 and 0.

    On the other hand, I4_MODP(A,360) is between 0 and 360, always.

  Example:

        I         J     MOD  I4_MODP   I4_MODP Factorization

      107        50       7       7    107 =  2 *  50 + 7
      107       -50       7       7    107 = -2 * -50 + 7
     -107        50      -7      43   -107 = -3 *  50 + 43
     -107       -50      -7      43   -107 =  3 * -50 + 43

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 January 2007

  Author:

    John Burkardt

  Parameters:

    Input, int I, the number to be divided.

    Input, int J, the number that divides I.

    Output, int I4_MODP, the nonnegative remainder when I is
    divided by J.
*/
{
  int value;

  if ( j == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "I4_MODP - Fatal error!\n" );
    fprintf ( stderr, "  I4_MODP ( I, J ) called with J = %d\n", j );
    exit ( 1 );
  }

  value = i % j;

  if ( value < 0 )
  {
    value = value + abs ( j );
  }

  return value;
}
/******************************************************************************/

int i4_wrap ( int ival, int ilo, int ihi )

/******************************************************************************/
/*
  Purpose:

    I4_WRAP forces an I4 to lie between given limits by wrapping.

  Example:

    ILO = 4, IHI = 8

    I   Value

    -2     8
    -1     4
     0     5
     1     6
     2     7
     3     8
     4     4
     5     5
     6     6
     7     7
     8     8
     9     4
    10     5
    11     6
    12     7
    13     8
    14     4

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 December 2012

  Author:

    John Burkardt

  Parameters:

    Input, int IVAL, an integer value.

    Input, int ILO, IHI, the desired bounds for the integer value.

    Output, int I4_WRAP, a "wrapped" version of IVAL.
*/
{
  int jhi;
  int jlo;
  int value;
  int wide;

  if ( ilo < ihi )
  {
   jlo = ilo;
   jhi = ihi;
  }
  else
  {
    jlo = ihi;
    jhi = ilo;
  }

  wide = jhi + 1 - jlo;

  if ( wide == 1 )
  {
    value = jlo;
  }
  else
  {
    value = jlo + i4_modp ( ival - jlo, wide );
  }

  return value;
}
/******************************************************************************/

void i4mat_transpose_print ( int m, int n, int a[], char *title )

/******************************************************************************/
/*
  Purpose:

    I4MAT_TRANSPOSE_PRINT prints an I4MAT, transposed.

  Discussion:

    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 January 2005

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, int A[M*N], the M by N matrix.

    Input, char *TITLE, a title.
*/
{
  i4mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
/******************************************************************************/

void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo,
  int ihi, int jhi, char *title )

/******************************************************************************/
/*
  Purpose:

    I4MAT_TRANSPOSE_PRINT_SOME prints some of an I4MAT, transposed.

  Discussion:

    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows of the matrix.
    M must be positive.

    Input, int N, the number of columns of the matrix.
    N must be positive.

    Input, int A[M*N], the matrix.

    Input, int ILO, JLO, IHI, JHI, designate the first row and
    column, and the last row and column to be printed.

    Input, char *TITLE, a title.
*/
{
# define INCX 10

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );

  if ( m <= 0 || n <= 0 )
  {
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  (None)\n" );
    return;
  }
/*
  Print the columns of the matrix, in strips of INCX.
*/
  for ( i2lo = ilo; i2lo <= ihi; i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;
    i2hi = i4_min ( i2hi, m );
    i2hi = i4_min ( i2hi, ihi );

    fprintf ( stdout, "\n" );
/*
  For each row I in the current range...

  Write the header.
*/
    fprintf ( stdout, "  Row: " );
    for ( i = i2lo; i <= i2hi; i++ )
    {
      fprintf ( stdout, "%6d  ", i - 1 );
    }
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  Col\n" );
    fprintf ( stdout, "\n" );
/*
  Determine the range of the rows in this strip.
*/
    j2lo = i4_max ( jlo, 1 );
    j2hi = i4_min ( jhi, n );

    for ( j = j2lo; j <= j2hi; j++ )
    {
/*
  Print out (up to INCX) entries in column J, that lie in the current strip.
*/
      fprintf ( stdout, "%5d: ", j - 1 );
      for ( i = i2lo; i <= i2hi; i++ )
      {
        fprintf ( stdout, "%6d  ", a[i-1+(j-1)*m] );
      }
      fprintf ( stdout, "\n" );
    }
  }

  return;
# undef INCX
}
/******************************************************************************/

void i4mat_write ( char *output_filename, int m, int n, int table[] )

/******************************************************************************/
/*
  Purpose:

    I4MAT_WRITE writes an I4MAT file.

  Discussion:

    An I4MAT is an array of I4's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 June 2009

  Author:

    John Burkardt

  Parameters:

    Input, char *OUTPUT_FILENAME, the output filename.

    Input, int M, the spatial dimension.

    Input, int N, the number of points.

    Input, int TABLE[M*N], the data.
*/
{
  int i;
  int j;
  FILE *output;
/*
  Open the file.
*/
  output = fopen ( output_filename, "wt" );

  if ( !output )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "I4MAT_WRITE - Fatal error!\n" );
    fprintf ( stderr, "  Could not open the output file.\n" );
    exit ( 1 );
  }
/*
  Write the data.
*/
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      fprintf ( output, "  %d", table[i+j*m] );
    }
    fprintf ( output, "\n" );
  }
/*
  Close the file.
*/
  fclose ( output );

  return;
}
/******************************************************************************/

void i4vec_print ( int n, int a[], char *title )

/******************************************************************************/
/*
  Purpose:

    I4VEC_PRINT prints an I4VEC.

  Discussion:

    An I4VEC is a vector of I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 November 2003

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of components of the vector.

    Input, int A[N], the vector to be printed.

    Input, char *TITLE, a title.
*/
{
  int i;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );
  fprintf ( stdout, "\n" );

  for ( i = 0; i < n; i++ )
  {
    fprintf ( stdout, "  %6d: %8d\n", i, a[i] );
  }
  return;
}
/******************************************************************************/

void interp ( char *code, int element_order, double r, double s, 
  double ubase[], double *u, double *dudr, double *duds )

/******************************************************************************/
/*
  Purpose: 

    INTERP interpolates a quantity in an element from basis node values.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    22 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, char *CODE, identifies the element.
    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 
    'T3', 'T6' and 'T10'.

    Input, int ELEMENT_ORDER, order of the element.

    Input, double R, S, the reference coordinates of a point.

    Input, double UBASE[ELEMENT_ORDER], the value of the quantity 
    at the basis nodes.

    Output, double *U, *DUDR, *DUDS, the interpolated value of the
    quantity and its derivatives at the point (R,S).
*/
{
  double *dtdr;
  double *dtds;
  int i;
  double *t;

  dtdr = ( double * ) malloc ( element_order * sizeof ( double ) );
  dtds = ( double * ) malloc ( element_order * sizeof ( double ) );
  t = ( double * ) malloc ( element_order * sizeof ( double ) );

  shape ( code, r, s, t, dtdr, dtds );
 
  *u = 0.0;
  for ( i = 0; i < element_order; i++ )
  {
    *u = *u + ubase[i] * t[i];
  }

  *dudr = 0.0;
  for ( i = 0; i < element_order; i++ )
  {
    *dudr = *dudr + ubase[i] * dtdr[i];
  }

  *duds = 0.0;
  for ( i = 0; i < element_order; i++ )
  {
    *duds = *duds + ubase[i] * dtds[i];
  }

  free ( dtdr );
  free ( dtds );
  free ( t );
 
  return;
}
/******************************************************************************/

void interp_test ( char *code )

/******************************************************************************/
/*
  Purpose:

    INTERP_TEST tests the interpolation property of an element.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    22 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, char *CODE, identifies the element.
    Legal values include "Q4", "Q8", "Q9", "Q12", "Q16", "QL", 
    "T3", "T4", "T6" and "T10".
*/
{
  double area;
  double dudr;
  double dudr_exact;
  double duds;
  double duds_exact;
  int element_order;
  int i;
  int node;
  double *node_r;
  double *node_s;
  double *node_u;
  double r;
  double r_factor;
  int *rexp;
  double s;
  double s_factor;
  int *sexp;
  int seed;
  int test;
  int test_num = 5;
  double u;
  double u_exact;

  if ( s_eqi ( code, "T4" ) )
  {
    printf ( "\n" );
    printf ( "INTERP_TEST - Warning!\n" );
    printf ( "  Skipping test for element \"T4\".\n" );
    return;
  }

  printf ( "\n" );
  printf ( "INTERP_TEST for element \"%s\".\n", code );

  element_order = order_code ( code );

  printf ( "  Number of nodes = %d\n", element_order );

  node_r = ( double * ) malloc ( element_order * sizeof ( double ) );
  node_s = ( double * ) malloc ( element_order * sizeof ( double ) );
  node_u = ( double * ) malloc ( element_order * sizeof ( double ) );
  rexp = ( int * ) malloc ( element_order * sizeof ( int ) );
  sexp = ( int * ) malloc ( element_order * sizeof ( int ) );
/*
  Get the coordinates of the reference nodes.
*/
  node_reference ( code, node_r, node_s, &area );
/*
  Get the monomial exponents for which the element is exact.
*/
  poly ( code, rexp, sexp );

  seed = 123456789;

  for ( i = 0; i < element_order; i++ )
  {
    printf ( "\n" );
    printf ( "  Interpolate R^%d * S^%d\n", rexp[i], sexp[i] );
    printf ( "\n" );
/*
  Evaluate R^REXP(I) * S^SEXP(I) at the nodes.  This is our data.
*/
    for ( node = 0; node < element_order; node++ )
    {
      r = node_r[node];
      s = node_s[node];
      if ( rexp[i] == 0 )
      {
        r_factor = 1.0;
      }
      else
      {
        r_factor = pow ( r, rexp[i] );
      }
      if ( sexp[i] == 0 )
      {
        s_factor = 1.0;
      }
      else
      {
        s_factor = pow ( s, sexp[i] );
      }
      node_u[node] = r_factor * s_factor;
      printf ( "  (R,S,U):       %12g  %12g  %12g\n", r, s, node_u[node] );
    }
/*
  Now pick random points in the element, and compute the interpolated
  value of R^REXP(*) * S^SEXP(I) there.  Mathematically, these
  values should be exact.
*/
    for ( test = 1; test <= test_num; test++ )
    {
      reference_sample ( code, &seed, &r, &s );

      printf ( "\n" );
      printf ( "  (R,S):            %12g  %12g\n", r, s );

      u_exact = r8_power ( r, rexp[i] ) * r8_power ( s, sexp[i] );

      dudr_exact = ( double ) ( rexp[i] ) 
        * r8_power ( r, rexp[i] - 1 ) * r8_power ( s, sexp[i] );

      duds_exact = r8_power ( r, rexp[i] ) * ( double ) ( sexp[i] ) 
        * r8_power ( s, sexp[i] - 1 );

      interp ( code, element_order, r, s, node_u, &u, &dudr, &duds );

      printf ( "  (U ,U* ,Error):   %12g  %12g  %12g\n", 
        u_exact, u, fabs ( u_exact - u ) );

      printf ( "  (Ur,Ur*,Error):   %12g  %12g  %12g\n",
        dudr_exact, dudr, fabs ( dudr_exact - dudr ) );

      printf ( "  (Us,Us*,Error):   %12g  %12g  %12g\n", 
        duds_exact, duds, fabs ( duds_exact - duds ) );
    }
  }

  free ( node_r );
  free ( node_s );
  free ( node_u ); 
  free ( rexp );
  free ( sexp );

  return;
}
/******************************************************************************/

void legendre_compute_dr ( int order, double xtab[], double weight[] )

/******************************************************************************/
/*
  Purpose:
  
    LEGENDRE_COMPUTE_DR: Gauss-Legendre quadrature by Davis-Rabinowitz method.
  
  Discussion:
    
    The integral:
  
      Integral ( -1 <= X <= 1 ) F(X) dX
  
    The quadrature rule:
  
      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    28 August 2007
  
  Author:
  
    Original FORTRAN77 version by Philip Davis, Philip Rabinowitz.
    C version by John Burkardt.
  
  Reference:
  
    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.
  
  Parameters:
  
    Input, int ORDER, the order.
    ORDER must be greater than 0.
  
    Output, double XTAB[ORDER], the abscissas.
  
    Output, double WEIGHT[ORDER], the weights.
*/
{
  double d1;
  double d2pn;
  double d3pn;
  double d4pn;
  double dp;
  double dpn;
  double e1;
  double fx;
  double h;
  int i;
  int iback;
  int k;
  int m;
  int mp1mi;
  int ncopy;
  int nmove;
  double p;
  double pi = 3.141592653589793;
  double pk;
  double pkm1;
  double pkp1;
  double t;
  double u;
  double v;
  double x0;
  double xtemp;

  if ( order < 1 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "LEGENDRE_COMPUTE_DR - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of ORDER = %d\n", order );
    exit ( 1 );
  }

  e1 = ( double ) ( order * ( order + 1 ) );

  m = ( order + 1 ) / 2;

  for ( i = 1; i <= m; i++ )
  {
    mp1mi = m + 1 - i;

    t = ( double ) ( 4 * i - 1 ) * pi / ( double ) ( 4 * order + 2 );

    x0 = cos ( t ) * ( 1.0 - ( 1.0 - 1.0 / ( double ) ( order ) )
      / ( double ) ( 8 * order * order ) );

    pkm1 = 1.0;
    pk = x0;