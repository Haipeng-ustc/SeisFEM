
int mesh_q16_element(int nelemx, int nelemy, int *element_node)

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


*/
{
  int base;
  int element;
  int element_order = 16;
  int i;
  int j;

  element = 0;

  for (j = 1; j <= nelemy; j++)
  {
    for (i = 1; i <= nelemx; i++)
    {
      base = (j - 1) * 3 * (3 * nelemx + 1) + 3 * i - 2;

      element_node[0 + element * element_order] = base;
      element_node[1 + element * element_order] = base + 1;
      element_node[2 + element * element_order] = base + 2;
      element_node[3 + element * element_order] = base + 3;
      element_node[4 + element * element_order] = base + (3 * nelemx + 1);
      element_node[5 + element * element_order] = base + (3 * nelemx + 1) + 1;
      element_node[6 + element * element_order] = base + (3 * nelemx + 1) + 2;
      element_node[7 + element * element_order] = base + (3 * nelemx + 1) + 3;
      element_node[8 + element * element_order] = base + 2 * (3 * nelemx + 1);
      element_node[9 + element * element_order] = base + 2 * (3 * nelemx + 1) + 1;
      element_node[10 + element * element_order] = base + 2 * (3 * nelemx + 1) + 2;
      element_node[11 + element * element_order] = base + 2 * (3 * nelemx + 1) + 3;
      element_node[12 + element * element_order] = base + 3 * (3 * nelemx + 1);
      element_node[13 + element * element_order] = base + 3 * (3 * nelemx + 1) + 1;
      element_node[14 + element * element_order] = base + 3 * (3 * nelemx + 1) + 2;
      element_node[15 + element * element_order] = base + 3 * (3 * nelemx + 1) + 3;

      element = element + 1;
    }
  }

}
