#include "dgemm.h"

#include <chrono>
#include <array>
#include <omp.h>
#include <iostream>
#include <string>
#include <sys/time.h>

extern "C"
{
#include "variorum.h"
}

/*
A, B    = input matrix (N x N)
C       = output matrix (N x N)
alpha   = scaling factor for matrix a
beta    = scaling factor for matrix c
*/

#ifndef N
#define N 1024
#elif N < 1024
#undef N
#define N 1024
#endif

#ifndef ALPHA
#define ALPHA 1.0
#endif

static double A[N][N];
static double B[N][N];
static double C[N][N];

/*
- still need some research how to set the appropriate N.
- limiting core frequency with variorum seems like doesn't work.
*/

void init_array()
{
#pragma omp parallel for
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            A[i][j] = ((double)i * j) / N;
            B[i][j] = ((double)i * j) / N;
        }
    }
}

void d_main()
{
    int t;
    std::cout << "Array size: " << N << std::endl;

#ifdef _OPENMP
    printf(HLINE);
#pragma omp parallel
    {
#pragma omp master
        {
            omp_display_affinity("Thread Affinity: %0.3L %.8n %.15{thread_affinity} %.12H");
            t = omp_get_num_threads();
            printf("Number of Threads requested = %i\n", t);
        }
    }
#endif

#ifdef _OPENMP
    t = 0;
#pragma omp parallel
#pragma omp critical
    t++;
    printf("Number of Threads counted = %i\n", t);
#endif

    init_array();

#pragma omp parallel for
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            C[i][j] = 0;

            for (int k = 0; k < N; k++)
            {
                C[i][j] += ALPHA * A[i][k] * B[k][j];
            }
        }
    }
}

double mysecond()
{
	struct timeval tp;
	struct timezone tzp;

	gettimeofday(&tp, &tzp);
	return ((double)tp.tv_sec + (double)tp.tv_usec * 1.e-6);
}

int main(int argc, char *argv[])
{
    if (argc < 4)
    {
        std::cerr << "Need argument : <core/socket/node> <id> <result file path>" << std::endl;
        return 1;
    }

    std::string level = argv[1];
    std::string filename = "result/DGEMM Benchmark on " + level + " level" + ".csv";

    // Create csv file for the result
    FILE *file = fopen(argv[3], "a");
    if (file == NULL)
    {
        std::cerr << "Error opening file." << std::endl;
        return -1;
    }

    // go to the end of file
    fseek(file, 0, SEEK_END);

    // print header if not printed yet
    if (ftell(file) == 0)
    {
        fprintf(file, "%sId, Time(uS)\n", level.c_str());
        fflush(file);
    }

    // run DGEMM benchmark
    /*     int ret;
    ret = variorum_cap_each_core_frequency_limit(1200);
    if (ret != 0)
    {
        printf("Print power failed!\n");
    } */
    int ret;

    ret = variorum_poll_power(stdout);
    double time = mysecond();

    d_main();

    time = mysecond()-time;
    ret = variorum_poll_power(stdout);
    //variorum_poll_power();
    /*         int ret = variorum_print_verbose_frequency();
            if (ret != 0)
            {
                printf("Print power failed!\n");
            } */

    // std::cout << res.count() << std::endl;

    fprintf(file, "%s,%f\n", argv[2], time);
    fflush(file);

    fclose(file);
    return 0;
}