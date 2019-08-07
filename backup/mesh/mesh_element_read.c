void mesh_element_read(char *element_name, int element_num, int element_order, int *element_node)
{
    int i, j;
    FILE *fp_element_node;
    if ((fp_element_node = fopen( element_name, "r")) == NULL)
        printf("\n element file cannot open\n");

    for (i = 0; i < element_num * element_order; i = i + element_order)
    {
        for (j = 0; j < element_order; j++)
        {
            fscanf(fp_element_node, "%d ", &element_node[i + j]);
        }
        fscanf(fp_element_node, "\n");
    }
    fclose(fp_element_node);

}
