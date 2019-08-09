
char *solver_type(int i)

/******************************************************************************/
/*
  Purpose:

    SOLVER_TYPE returns the code for liner system solver.

  List:

   SOLVER_TYPE List:
    I  SOLVER_TYPE   Definition
    -  ------------   ----------
    1  pardiso        require license and only valid for username haipeng;
    2  mgmres         Generalized Minimum Residual (GMRES) algorithm, CSR format;
    3  masslump       masslump.

*/
{
  char *value;

  if (i == 1)
  {
    value = (char *)malloc(8 * sizeof(char));
    strcpy(value, "pardiso");
  }
  else if (i == 2)
  {
    value = (char *)malloc(8 * sizeof(char));
    strcpy(value, "mgmres");
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
