
void stiffness_sparse_q16(int node_num, int element_num, int element_order, int *element_node, double **node_xy, double *vp, int *stiffi, int *stiffj, double *stiffness)

/* stiffness_sparse_q16 computes the stiffness matrix, store the matrix in the coo format (i,j,value), using 16-node rectangular, 36 ponits quadrature rule.

  Reference Element Q16:

    |
    1 13--14--15--16
    |  |   :   :   |
    |  |   :   :   |
    |  9..10..11..12
    S  |   :   :   |
    |  |   :   :   |
    |  5...6...7...8
    |  |   :   :   |
    |  |   :   :   |  
    0  1---2---3---4
    |
    +--0-----R-----1-->

    ELEMENT_NODE = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16		

*/
{

  int quad_num = 36;
  int i, j, iq, jq, ip, jp, element, quad, coo_index;
  int p1, p2, p3, p4;
  double r, s;
  double rtab[36], stab[36], weight[36];
  double w[16], dwdr[16], dwds[16];
  double phi[16], dphidx[16], dphidy[16];
  double x1, x2, x3, x4;
  double y1, y2, y3, y4;
  double area, det;
  double dxdr, dxds, dydr, dyds;
  double drdx, drdy, dsdx, dsdy;
  double a, b, c, d, e, f;

  a = 0.033765242898424;
  b = 0.169395306766868;
  c = 0.380690406958402;
  d = 0.619309593041599;
  e = 0.830604693233132;
  f = 0.966234757101576;

  rtab[0] = a;
  rtab[1] = a;
  rtab[2] = a;
  rtab[3] = a;
  rtab[4] = a;
  rtab[5] = a;
  rtab[6] = b;
  rtab[7] = b;
  rtab[8] = b;
  rtab[9] = b;
  rtab[10] = b;
  rtab[11] = b;
  rtab[12] = c;
  rtab[13] = c;
  rtab[14] = c;
  rtab[15] = c;
  rtab[16] = c;
  rtab[17] = c;
  rtab[18] = d;
  rtab[19] = d;
  rtab[20] = d;
  rtab[21] = d;
  rtab[22] = d;
  rtab[23] = d;
  rtab[24] = e;
  rtab[25] = e;
  rtab[26] = e;
  rtab[27] = e;
  rtab[28] = e;
  rtab[29] = e;
  rtab[30] = f;
  rtab[31] = f;
  rtab[32] = f;
  rtab[33] = f;
  rtab[34] = f;
  rtab[35] = f;

  stab[0] = a;
  stab[1] = b;
  stab[2] = c;
  stab[3] = d;
  stab[4] = e;
  stab[5] = f;
  stab[6] = a;
  stab[7] = b;
  stab[8] = c;
  stab[9] = d;
  stab[10] = e;
  stab[11] = f;
  stab[12] = a;
  stab[13] = b;
  stab[14] = c;
  stab[15] = d;
  stab[16] = e;
  stab[17] = f;
  stab[18] = a;
  stab[19] = b;
  stab[20] = c;
  stab[21] = d;
  stab[22] = e;
  stab[23] = f;
  stab[24] = a;
  stab[25] = b;
  stab[26] = c;
  stab[27] = d;
  stab[28] = e;
  stab[29] = f;
  stab[30] = a;
  stab[31] = b;
  stab[28] = c;
  stab[29] = d;
  stab[30] = e;
  stab[31] = f;

  weight[0] = 0.007338020422245;
  weight[1] = 0.015451823343096;
  weight[2] = 0.020041279329452;
  weight[3] = 0.020041279329452;
  weight[4] = 0.015451823343096;
  weight[5] = 0.007338020422245;
  weight[6] = 0.015451823343096;
  weight[7] = 0.032537228147042;
  weight[8] = 0.042201341771897;
  weight[9] = 0.042201341771897;
  weight[10] = 0.032537228147042;
  weight[11] = 0.015451823343096;
  weight[12] = 0.020041279329452;
  weight[13] = 0.042201341771897;
  weight[14] = 0.054735862541824;
  weight[15] = 0.054735862541824;
  weight[16] = 0.042201341771897;
  weight[17] = 0.020041279329452;
  weight[18] = 0.020041279329452;
  weight[19] = 0.042201341771897;
  weight[20] = 0.054735862541824;
  weight[21] = 0.054735862541824;
  weight[22] = 0.042201341771897;
  weight[23] = 0.020041279329452;
  weight[24] = 0.015451823343096;
  weight[25] = 0.032537228147042;
  weight[26] = 0.042201341771897;
  weight[27] = 0.042201341771897;
  weight[28] = 0.032537228147042;
  weight[29] = 0.015451823343096;
  weight[30] = 0.007338020422245;
  weight[31] = 0.015451823343096;
  weight[32] = 0.020041279329452;
  weight[33] = 0.020041279329452;
  weight[34] = 0.015451823343096;
  weight[35] = 0.007338020422245;

  coo_index = 0; // for coo format sparse matrix index
  for (i = 0; i < element_num * element_order * element_order; i++)
  {
    stiffi[i] = 0;
    stiffj[i] = 0;
    stiffness[i] = 0.0;
  }

  for (element = 0; element < element_num * element_order; element = element + element_order)
  {
    p1 = element_node[0 + element] - 1;
    p2 = element_node[3 + element] - 1;
    p3 = element_node[12 + element] - 1;
    p4 = element_node[15 + element] - 1;
    x1 = node_xy[0][p1];
    y1 = node_xy[1][p1];
    x2 = node_xy[0][p2];
    y2 = node_xy[1][p2];
    x3 = node_xy[0][p3];
    y3 = node_xy[1][p3];
    x4 = node_xy[0][p4];
    y4 = node_xy[1][p4];

    area = 0.5 * fabs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) + 0.5 * fabs(x1 * (y4 - y3) + x4 * (y3 - y1) + x3 * (y1 - y4));

    if (area == 0.0)
    {
      printf("STIFFNESS_SPARSE_q16 - Fatal error!\n");
      printf("Zero area for element: %d\n", element);
      exit(1);
    }

    //  For each quadrature point in the element...

    for (quad = 0; quad < quad_num; quad++)
    {
      r = rtab[quad];
      s = stab[quad];

      dxdr = -(1 - s) * x1 + (1 - s) * x2 + s * x3 - s * x4;
      dxds = -(1 - r) * x1 - r * x2 + r * x3 + (1 - r) * x4;
      dydr = -(1 - s) * y1 + (1 - s) * y2 + s * y3 - s * y4;
      dyds = -(1 - r) * y1 - r * y2 + r * y3 + (1 - r) * y4;

      det = dxdr * dyds - dxds * dydr;
      drdx = dyds / det;
      drdy = -dxds / det;
      dsdx = -dydr / det;
      dsdy = dxdr / det;

      shape_q16(r, s, w, dwdr, dwds);

      phi[0] = w[0];
      phi[1] = w[1];
      phi[2] = w[2];
      phi[3] = w[3];
      phi[4] = w[4];
      phi[5] = w[5];
      phi[6] = w[6];
      phi[7] = w[7];
      phi[8] = w[8];
      phi[9] = w[9];
      phi[10] = w[10];
      phi[11] = w[11];
      phi[12] = w[12];
      phi[13] = w[13];
      phi[14] = w[14];
      phi[15] = w[15];

      dphidx[0] = dwdr[0] * drdx + dwds[0] * dsdx;
      dphidx[1] = dwdr[1] * drdx + dwds[1] * dsdx;
      dphidx[2] = dwdr[2] * drdx + dwds[2] * dsdx;
      dphidx[3] = dwdr[3] * drdx + dwds[3] * dsdx;
      dphidx[4] = dwdr[4] * drdx + dwds[4] * dsdx;
      dphidx[5] = dwdr[5] * drdx + dwds[5] * dsdx;
      dphidx[6] = dwdr[6] * drdx + dwds[6] * dsdx;
      dphidx[7] = dwdr[7] * drdx + dwds[7] * dsdx;
      dphidx[8] = dwdr[8] * drdx + dwds[8] * dsdx;
      dphidx[9] = dwdr[9] * drdx + dwds[9] * dsdx;
      dphidx[10] = dwdr[10] * drdx + dwds[10] * dsdx;
      dphidx[11] = dwdr[11] * drdx + dwds[11] * dsdx;
      dphidx[12] = dwdr[12] * drdx + dwds[12] * dsdx;
      dphidx[13] = dwdr[13] * drdx + dwds[13] * dsdx;
      dphidx[14] = dwdr[14] * drdx + dwds[14] * dsdx;
      dphidx[15] = dwdr[15] * drdx + dwds[15] * dsdx;

      dphidy[0] = dwdr[0] * drdy + dwds[0] * dsdy;
      dphidy[1] = dwdr[1] * drdy + dwds[1] * dsdy;
      dphidy[2] = dwdr[2] * drdy + dwds[2] * dsdy;
      dphidy[3] = dwdr[3] * drdy + dwds[3] * dsdy;
      dphidy[4] = dwdr[4] * drdy + dwds[4] * dsdy;
      dphidy[5] = dwdr[5] * drdy + dwds[5] * dsdy;
      dphidy[6] = dwdr[6] * drdy + dwds[6] * dsdy;
      dphidy[7] = dwdr[7] * drdy + dwds[7] * dsdy;
      dphidy[8] = dwdr[8] * drdy + dwds[8] * dsdy;
      dphidy[9] = dwdr[9] * drdy + dwds[9] * dsdy;
      dphidy[10] = dwdr[10] * drdy + dwds[10] * dsdy;
      dphidy[11] = dwdr[11] * drdy + dwds[11] * dsdy;
      dphidy[12] = dwdr[12] * drdy + dwds[12] * dsdy;
      dphidy[13] = dwdr[13] * drdy + dwds[13] * dsdy;
      dphidy[14] = dwdr[14] * drdy + dwds[14] * dsdy;
      dphidy[15] = dwdr[15] * drdy + dwds[15] * dsdy;

      for (iq = 0; iq < element_order; iq++)
      {
        ip = element_node[iq + element] - 1; // c array from 0

        for (jq = 0; jq < element_order; jq++)
        {
          jp = element_node[jq + element] - 1; // c array from 0

          coo_index = element * element_order + iq * element_order + jq;
          stiffi[coo_index] = ip;
          stiffj[coo_index] = jp;
          stiffness[coo_index] = stiffness[coo_index] + area * vp[ip] * weight[quad] * (dphidx[iq] * dphidx[jq] + dphidy[iq] * dphidy[jq]);
        }
      }
    }
  }
}
