
char *element_type ( int i )

/******************************************************************************/
/*
  Purpose:

    ELEMENT_TYPE returns the code for each element.

  List:

    I  ELEMENT_TYPE   Definition
    -  ------------   ----------
    1  T3             3 node linear triangle;
    2  T6             6 node quadratic triangle;
    3  T10            10 node cubic triangle.
    4  Q4             4 node linear Lagrange/serendipity quadrilateral;
    5  Q9             9 node quadratic Lagrange quadrilateral;
    6  Q16            16 node cubic Lagrange quadrilateral;

*/
{
  char *value;

  if ( i == 1 )
  {
    value = ( char * ) malloc ( 3 * sizeof ( char ) );
    strcpy ( value, "T3" );
  }
  else if ( i == 2 )
  {
    value = ( char * ) malloc ( 3 * sizeof ( char ) );
    strcpy ( value, "T6" );
  }
  else if ( i == 3 )
  {
    value = ( char * ) malloc ( 4 * sizeof ( char ) );
    strcpy ( value, "T10" );
  }
  else if ( i == 4 )
  {
    value = ( char * ) malloc ( 3 * sizeof ( char ) );
    strcpy ( value, "Q4" );
  }
  else if ( i == 5 )
  {
    value = ( char * ) malloc ( 3 * sizeof ( char ) );
    strcpy ( value, "Q9" );
  }
  else if ( i == 6 )
  {
    value = ( char * ) malloc ( 4 * sizeof ( char ) );
    strcpy ( value, "Q16" );
  }
  else
  {
    value = ( char * ) malloc ( 4 * sizeof ( char ) );
    strcpy ( value, "???" );
    printf("ELEMENT_TYPE - Fatal error!\n");
    exit(1);
  }

  return value;
}
