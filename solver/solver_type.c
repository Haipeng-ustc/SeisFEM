
char *solver_type(int i)

/******************************************************************************/
/*
  Purpose:

    SOLVER_TYPE returns the code for liner system solver.

  List:

    I  SOLVER_TYPE   Definition
    -  ------------   ----------
    1  mgmres         Generalized Minimum Residual (GMRES) algorithm,
                      using compressed row sparse matrix format
    2  superlu        LU decomposition to solver linear system;
    3  masslump       masslump.

*/
{
  char *value;

  if (i == 1)
  {
    value = (char *)malloc(6 * sizeof(char));
    strcpy(value, "mgmres");
  }
  else if (i == 2)
  {
    value = (char *)malloc(7 * sizeof(char));
    strcpy(value, "superlu");
  }
  else if (i == 3)
  {
    value = (char *)malloc(8 * sizeof(char));
    strcpy(value, "masslump");
  }
  else
  {
    value = (char *)malloc(4 * sizeof(char));
    strcpy(value, "???");
    printf("SOLVER_TYPE - Fatal error!\n");
    exit(1);
  }

  return value;
}
