#include "ranker_aux.hpp"

#include <chrono>
#include <array>
#include <omp.h>
#include <iostream>
#include <string>
#include <sys/time.h>
#include <thread>
#include <atomic>

/*
A, B    = input matrix (N x N)
C       = output matrix (N x N)
alpha   = scaling factor for matrix a
beta    = scaling factor for matrix c
*/
#define HLINE "-------------------------------------------------------------\n"

#ifndef N
#define N 1024
#endif

#ifndef ALPHA
#define ALPHA 1.5
#endif

#ifndef NTIMES
#define NTIMES 11
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
            C[i][j] = 0.0;
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

void d_main(const std::string &node, const std::string &socket, std::vector<double> &time, std::vector<double> &energy, std::vector<double> &power)
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

    for (int iter = 0; iter < NTIMES; iter++)
    {
        record_power("dgemm_" + node, "w");

        double clock = mysecond();

#pragma omp parallel for
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                for (int k = 0; k < N; k++)
                {
                    C[i][j] += ALPHA * A[i][k] * B[k][j];
                }
            }
        }

        clock = mysecond() - clock;
        time.push_back(clock);

        record_power("dgemm_" + node, "a");

        std::pair<double, double> parse = parse_poll_power_pkg("dgemm_" + node, socket);
        energy.push_back(parse.first);
        power.push_back(parse.second);
    }
    return;
}

int main(int argc, char *argv[])
{
    if (argc < 7)
    {
        std::cerr << "Need argument : <Node Id> <Socket Id> <Core Id> <Idle Power> <Power limit> <result file path>" << std::endl;
        return 1;
    }

    std::string node = argv[1];
    std::string socket = argv[2];
    std::string core = argv[3];
    std::string power_limit = argv[5];

    // Create csv file for the result
    FILE *file = fopen(argv[6], "a");
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
        fprintf(file, "NodeId,SocketId,CPUId,NumberOfThreads,PowerLimit(W),IdleCPUPower(W),AvgTime(S),SDTime,AvgEnergyUsage(J),SDEnergy,AvgPowerConsumption(W),SDPowerConsumption,TotalFLOP,AvgFLOPS,SDFLOPS,AvguJ/FLOP,SDuJ/FLOP\n");
        fflush(file);
    }

    // MAIN
    std::vector<double> time, energy, power_cons;
    int n_threads;

#pragma omp parallel
    {
#pragma omp master
        {
            n_threads = omp_get_num_threads();
        }
    }
    d_main(node, socket, time, energy, power_cons);

    // Remove first element
    time.erase(time.begin());
    energy.erase(energy.begin());
    power_cons.erase(power_cons.begin());

    // Time
    double avg_time, sd_time;
    standard_deviation(time, avg_time, sd_time);

    // Energy
    double avg_energy, sd_energy;
    standard_deviation(energy, avg_energy, sd_energy);

    // Power Consumption
    double avg_power, sd_power;
    standard_deviation(power_cons, avg_power, sd_power);

    double idle_power = stod(std::string(argv[4]));

    // Flop
    u_int64_t n64 = static_cast<u_int64_t>(N);
    u_int64_t flop = 3 * n64 * n64 * n64;

    std::vector<double> flops, jpflop;
    for (size_t i = 0; i < time.size(); i++)
    {
        flops.push_back(flop / time[i]);
        jpflop.push_back(energy[i] * 1000000 / flop);
    }

    // FLOPS
    double avg_flops, sd_flops;
    standard_deviation(flops, avg_flops, sd_flops);

    // FLOPS/W
    double avg_jpflop, sd_jpflop;
    standard_deviation(jpflop, avg_jpflop, sd_jpflop);

    //"NodeId,SocketId,ThreadId,NumberOfThreads,PowerLimit(W),IdleCPUPower(W),AvgTime(S),SDTime,AvgEnergyUsage(J),SDEnergy,AvgPowerConsumption(W),SDPowerConsumption,TotalFLOP,AvgFLOPS,SDFLOPS,AvguJ/FLOP,SDuJ/FLOP\n"
    fprintf(file, "%s,%s,%s,%d,%s,%f,%f,%f,%f,%f,%f,%f,%lu,%f,%f,%f,%f\n",
            node.c_str(), socket.c_str(), core.c_str(), n_threads, power_limit.c_str(), idle_power, avg_time, sd_time, avg_energy, sd_energy, avg_power, sd_power, flop, avg_flops, sd_flops, avg_jpflop, sd_jpflop);
    fflush(file);

    fclose(file);
    return 0;
}