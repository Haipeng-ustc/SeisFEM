#include <stdio.h>

int main()
{   
    FILE *fp;
    int i ,node_num;
    double data1, data2, data3, data4, data5, data6, data7;
    node_num =  104974;
    data1 = 2200.000000;
    data2 = 2000.000000;
    data3 = 1154.700000;
    data4 = 8800000000.000000;
    data5 = 2933338804.000000;
    data6 = 8800000000.000000;
    data7 = 2933330598.000000;
    fp = fopen("velocity_and_density.txt","w");
    for(i=0;i<node_num ;i++)
    {
        fprintf(fp, "%f %f  %f  %f  %f  %f  %f\n",data1,data2, data3, data4, data5, data6, data7);
    }
    fclose(fp);
 



}