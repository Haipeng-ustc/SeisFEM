printf("\n Write some arrays into files. \n");
if ((fp_element_node = fopen("element_node.txt", "w")) == NULL)
      printf("\n File: element_node.txt cannot open\n");
else
{
      for (i = 0; i < element_num * element_order; i = i + element_order)
      {
            for (order = 0; order < element_order; order++)
            {
                  fprintf(fp_element_node, "%d	", element_node[order + i]);
            }
            fprintf(fp_element_node, "\n");
      }
      printf(" Array: element_node has ben written into file element_node.txt\n");
}

if ((fp_node_xy = fopen("node_xy.txt", "w")) == NULL)
      printf("\n File: node_xy.txt cannot open\n");
else
{
      for (i = 0; i < node_num; i++)
      {
            fprintf(fp_node_xy, "%f, %f\n", node_xy[0][i], node_xy[1][i]);
      }
      printf(" Array: node_xy has been written into file node_xy.txt\n");
}

/* if (( fp_mass = fopen("mass_matrix.txt","w")) == NULL) printf("\n File: mass_matrix.txt cannot open\n");
  else  
  {	for( i = 0; i < nnz; i++) fprintf(fp_mass, "%d	%d	%f\n",mi[i] + 1, mj[i] + 1, mass[i]);   // 0-based index to 1-based index 
	printf(" Array: mi, mj and mass have been written into file mass_matrix.txt\n");	  
  }
*/
if ((fp_mass = fopen("mass_matrix_lump.txt", "w")) == NULL)
      printf("\n File: mass_matrix_lump.txt cannot open\n");
else
{
      for (i = 0; i < node_num; i++)
      {
            fprintf(fp_mass, "%f\n", mass[i]); // 0-based index to 1-based index
      }
      printf(" Array: mass have been written into file mass_matrix_lump.txt\n");
}

if ((fp_stiffness = fopen("stiffness_matrix.txt", "w")) == NULL)
      printf("\n File: stiffness_matrix.txt cannot open\n");
else
{
      for (i = 0; i < nnz; i++)
            fprintf(fp_stiffness, "%d	%d	%f\n", stiffi[i] + 1, stiffj[i] + 1, stiffness[i]); // 0-based index to 1-based index
      printf(" Array: stiffi, stiffj and stiffness have been written into file stiffness_matrix.txt\n");
}
if ((fp_model = fopen("acoustic_model.txt", "w")) == NULL)
      printf("\n File: acoustic_model.txt cannot open\n");
else
{
      for (i = 0; i < node_num; i++)
      {
            fprintf(fp_model, "%f	%f	%f\n", rho[i], vp[i], vs[i]);
      }
      printf(" Array: acoustic_model has been written into file acoustic_model.txt\n");
}

if ((fp_model2 = fopen("elastic_model.txt", "w")) == NULL)
      printf("\n File: elastic_model.txt cannot open\n");
else
{
      for (i = 0; i < node_num; i++)
      {
            fprintf(fp_model2, "%f	%f	%f	%f	%f\n", rho[i], c11[i], c13[i], c33[i], c44[i]);
      }
      printf(" Array: elastic_model has been written into file elastic_model.txt\n");
}

if ((fp_pml = fopen("pml.txt", "w")) == NULL)
      printf("\n File: pml.txt cannot open\n");
else
{
      for (i = 0; i < node_num; i++)
      {
            fprintf(fp_pml, "%f	%f	%f	%f\n", pml_dx[i], pml_dy[i], pml_dxx[i], pml_dyy[i]);
      }
      printf(" Array: pml has been written into file pml.txt\n");
}

if ((fp_mpml = fopen("mpml.txt", "w")) == NULL)
      printf("\n File: mpml.txt cannot open\n");
else
{
      for (i = 0; i < node_num; i++)
      {
            fprintf(fp_mpml, "%f	%f	%f	%f	%f	%f\n", mpml_dx[i], mpml_dy[i], mpml_dxx[i], mpml_dyy[i], mpml_dxx_pyx[i], mpml_dyy_pxy[i]);
      }
      printf(" Array: mpml has been written into file mpml.txt\n");
}
