int mesh_node_num_read(char *mesh_par_name)
{
    int element_num, node_num;
    FILE *fp_mesh_par;
    if ((fp_mesh_par = fopen(mesh_par_name, "r")) == NULL)
        printf("\n mesh_par file cannot open\n");

    fscanf(fp_mesh_par, "element_num = %d\n", &element_num);
    fscanf(fp_mesh_par, "node_num = %d\n", &node_num);
    fclose(fp_mesh_par);

    return node_num;
}
