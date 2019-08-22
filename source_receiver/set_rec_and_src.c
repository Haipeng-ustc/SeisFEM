void set_rec_and_src(int rec_num, double rec_x_first, double rec_y_first, double rec_x_last, double rec_y_last, double *rec_x, double *rec_y,
                      int src_num, double src_x_first, double src_y_first, double src_x_last, double src_y_last, double *src_x, double *src_y)
{
    int i;
    for (i = 0; i < rec_num; i++)
    {
        if (rec_num == 1)
        {
            rec_x[i] = rec_x_first;
            rec_y[i] = rec_y_first;
        }
        else
        {
            rec_x[i] = rec_x_first + i * (rec_x_last - rec_x_first) / (1.0 * (rec_num - 1));
            rec_y[i] = rec_y_first + i * (rec_y_last - rec_y_first) / (1.0 * (rec_num - 1));
        }
    }

    for (i = 0; i < src_num; i++)
    {
        if (src_num == 1)
        {
            src_x[i] = src_x_first;
            src_y[i] = src_y_first;
        }
        else
        {
            src_x[i] = src_x_first + i * (src_x_last - src_x_first) / (1.0 * (src_num - 1));
            src_y[i] = src_y_first + i * (src_y_last - src_y_first) / (1.0 * (src_num - 1));
        }
    }

}