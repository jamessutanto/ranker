#include "dgemm.h"
#include "variorum_parser.hpp"

#include <chrono>
#include <array>
#include <omp.h>
#include <iostream>
#include <string>
#include <sys/time.h>

/*
A, B    = input matrix (N x N)
C       = output matrix (N x N)
alpha   = scaling factor for matrix a
beta    = scaling factor for matrix c
*/

#ifndef N
#define N 1024
#endif

#ifndef ALPHA
#define ALPHA 1.0
#endif

static double A[N][N];
static double B[N][N];
static double C[N][N];

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

double mysecond()
{
    struct timeval tp;
    struct timezone tzp;

    gettimeofday(&tp, &tzp);
    return ((double)tp.tv_sec + (double)tp.tv_usec * 1.e-6);
}

std::pair<double, double> d_main(const std::string &node, const std::string &socket)
{
    int t;
    std::cout << "Array size: " << N << std::endl;
    std::array<double, 3> time, energy;

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

    for (int k = 0; k < 3; k++)
    {
        record_power("dgemm_" + node, "w");
        time[k] = mysecond();
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
        time[k] = mysecond() - time[k];
        record_power("dgemm_" + node, "a");
        energy[k] = parse_poll_power_pkg("dgemm_" + node, socket);
    }
    double avgtime = 0;
    double usage = 0;
    for (int k = 0; k < 3; k++)
    {
        avgtime += time[k];
        usage += energy[k];
    }
    avgtime = avgtime / 3;
    usage = usage / 3;

    return std::make_pair(avgtime, usage);
}

int main(int argc, char *argv[])
{
    if (argc < 5)
    {
        std::cerr << "Need argument : <Node Id> <Socket Id> <Core Id> <Idle Power> <result file path>" << std::endl;
        return 1;
    }

    std::string node = argv[1];
    std::string socket = argv[2];
    std::string core = argv[3];

    // Create csv file for the result
    FILE *file = fopen(argv[5], "a");
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
        fprintf(file, "NodeId, SocketId, ThreadId, StaticEnergyUsage(J), Time(S), DynamicEnergyUsage(J)\n");
        fflush(file);
    }

    // MAIN
    record_power("dgemm_" + node, "w");

    std::pair<double, double> result = d_main(node, socket);

    double idle_power = stod(std::string(argv[4])) * result.first;    

    fprintf(file, "%s,%s,%s,%f,%f,%f\n", node.c_str(), socket.c_str(), core.c_str(), idle_power, result.first, result.second - idle_power);
    fflush(file);

    fclose(file);
    return 0;
}