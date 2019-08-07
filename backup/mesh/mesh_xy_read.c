void mesh_xy_read(char *node_xy_name, int node_num, double **node_xy)
{
    int i;
    double x, y;
    FILE *fp_node_xy;
    if ((fp_node_xy = fopen( node_xy_name, "r")) == NULL)
        printf("\n node_xy file cannot open\n");
    for (i = 0; i < node_num; i++)
    {
        fscanf(fp_node_xy, "%lf %lf\n", &x, &y);
        node_xy[0][i] = x;
        node_xy[1][i] = y;
    }
    fclose(fp_node_xy);
}
