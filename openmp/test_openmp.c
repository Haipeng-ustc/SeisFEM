#include <stdio.h>
#include <omp.h>

int main()
{
    int tid;    //表示当前线程号
    int nthreads;   //表示所有的线程数
    omp_set_num_threads(3);     //设置启用3个线程

    int a[3]={1,2,3};
    int b[3]={1,2,3};
    int c[3];

    /*并行区域开始*/
    #pragma omp parallel private(tid,nthreads) shared(a,b,c)
    {
        #pragma omp for
        for(int i=0;i<3;i++)
        {
            tid=omp_get_thread_num();   //获取当前的进程号
            nthreads=omp_get_num_threads();     //获取总进程数
            c[i]=a[i]+b[i];
        }
    }
    /*并行区域结束*/

    return 0;
}