#include "ranker.h"
#include "ranker_aux.hpp"

#include <iostream>
#include <cstring>
#include <hwloc.h>
#include <string>
#include <vector>
#include <sstream>
#include <sched.h>
#include <sys/wait.h>
#include <mpi.h>
#include <fstream>
#include <filesystem>
#include <cmath>
#include <fcntl.h>
#include <chrono>

extern "C"
{
#include "variorum.h"
}

#include "sys-sage.hpp"

void run(const std::string &, const std::string &, const std::string &, const int);

void rank_dgemm_core(const std::string &result)
{
    int powerLimit = get_socket_power_limit();

    // Do not change the number
    run("dgemm", "core", result + "-1,0", powerLimit);

    run("dgemm", "core", result + "-0,7", powerLimit);

    run("dgemm", "core", result + "-0,5", powerLimit);

    run("dgemm", "core", result + "-0,4", powerLimit);

    run("dgemm", "core", result + "-0,3", powerLimit);

    run("dgemm", "core", result + "-0,2", powerLimit);
}

void rank_dgemm_socket(const std::string &result)
{
    int powerLimit = get_socket_power_limit();

    // Do not change the number
    run("dgemm", "socket", result + "-1,0", powerLimit);

/*     run("dgemm", "socket", result + "-0,7", powerLimit);

    run("dgemm", "socket", result + "-0,5", powerLimit);

    run("dgemm", "socket", result + "-0,4", powerLimit);

    run("dgemm", "socket", result + "-0,3", powerLimit);

    run("dgemm", "socket", result + "-0,2", powerLimit); */
}

void rank_dgemm_node(const std::string &result)
{
    int powerLimit = get_node_power_limit();

    // Do not change the number
    run("dgemm", "node", result + "-1,0", powerLimit);

    run("dgemm", "node", result + "-0,7", powerLimit);

    run("dgemm", "node", result + "-0,5", powerLimit);

    run("dgemm", "node", result + "-0,4", powerLimit);

    run("dgemm", "node", result + "-0,3", powerLimit);

    run("dgemm", "node", result + "-0,2", powerLimit);
}

void rank_stream_socket(const std::string &result)
{
    run("stream", "socket", result, -1);
}

void rank_stream_node(const std::string &result)
{
    run("stream", "node", result, -1);
}

void rank_srandom_socket(const std::string &result)
{
    run("srandom", "socket", result, -1);
}

void rank_srandom_node(const std::string &result)
{
    run("srandom", "node", result, -1);
}

void calculate_rank(const std::string &level)
{
    std::vector<std::string> power = {"1,0", "0,7", "0,5", "0,4", "0,3", "0,2"};
    for (std::string s : power)
    {
        std::string ep = "result/data/ranker_" + level + "_ep_" + s + "_rank.csv";
        std::string bt = "result/data/ranker_" + level + "_bt_" + s + "_rank.csv";
        std::string sp = "result/data/ranker_" + level + "_sp_" + s + "_rank.csv";
    }
}

int execute(const std::string &command, const std::string &arg1, const std::string &arg2,
            const std::string &arg3, const std::string &arg4, const std::string &arg5, const std::string &arg6);

void rank_core()
{
    if (execute("make", "-C", "./NPB3.4.2/NPB3.4-OMP/", "ep", "CLASS=A", "", "") != 0)
    {
        std::cerr << "Failed to compile NAS Parallel EP benchmark" << std::endl;
        return;
    }

    if (execute("make", "-C", "./NPB3.4.2/NPB3.4-OMP/", "bt", "CLASS=A", "", "") != 0)
    {
        std::cerr << "Failed to compile NAS Parallel EP benchmark" << std::endl;
        return;
    }

    if (execute("make", "-C", "./NPB3.4.2/NPB3.4-OMP/", "sp", "CLASS=A", "", "") != 0)
    {
        std::cerr << "Failed to compile NAS Parallel EP benchmark" << std::endl;
        return;
    }

    MPI_Init(NULL, NULL);

    std::filesystem::path cwd = std::filesystem::current_path();

    rank_dgemm_core((cwd / "result/data/ranker_core_dgemm").string());

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    calculate_rank("core");
}

void rank_socket()
{
    if (execute("make", "-C", "./NPB3.4.2/NPB3.4-OMP/", "ep", "CLASS=C", "", "") != 0)
    {
        std::cerr << "Failed to compile NAS Parallel EP benchmark" << std::endl;
        return;
    }

    if (execute("make", "-C", "./NPB3.4.2/NPB3.4-OMP/", "bt", "CLASS=C", "", "") != 0)
    {
        std::cerr << "Failed to compile NAS Parallel EP benchmark" << std::endl;
        return;
    }

    if (execute("make", "-C", "./NPB3.4.2/NPB3.4-OMP/", "sp", "CLASS=C", "", "") != 0)
    {
        std::cerr << "Failed to compile NAS Parallel EP benchmark" << std::endl;
        return;
    }
    MPI_Init(NULL, NULL);

    std::filesystem::path cwd = std::filesystem::current_path();

    rank_stream_socket(cwd / "result/data/ranker_socket_stream");
    rank_dgemm_socket(cwd / "result/data/ranker_socket_dgemm");

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    calculate_rank("socket");
}

void rank_node()
{
    if (execute("make", "-C", "./NPB3.4.2/NPB3.4-OMP/", "ep", "CLASS=D", "", "") != 0)
    {
        std::cerr << "Failed to compile NAS Parallel EP benchmark" << std::endl;
        return;
    }

    if (execute("make", "-C", "./NPB3.4.2/NPB3.4-OMP/", "bt", "CLASS=D", "", "") != 0)
    {
        std::cerr << "Failed to compile NAS Parallel EP benchmark" << std::endl;
        return;
    }

    if (execute("make", "-C", "./NPB3.4.2/NPB3.4-OMP/", "sp", "CLASS=D", "", "") != 0)
    {
        std::cerr << "Failed to compile NAS Parallel EP benchmark" << std::endl;
        return;
    }

    MPI_Init(NULL, NULL);

    std::filesystem::path cwd = std::filesystem::current_path();

    rank_stream_node(cwd / "result/data/ranker_node_stream");
    rank_dgemm_node(cwd / "result/data/ranker_node_dgemm");

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    calculate_rank("node");
}

int init_topo(const char *filename)
{
    int err;
    unsigned long flags = 0; // don't show anything special
    hwloc_topology_t topology;

    err = hwloc_topology_init(&topology);
    if (err)
    {
        std::cerr << "hwloc: Failed to initialize" << std::endl;
        return 1;
    }
    err = hwloc_topology_set_flags(topology, flags);
    if (err)
    {
        std::cerr << "hwloc: Failed to set flags" << std::endl;
        return 1;
    }
    err = hwloc_topology_load(topology);
    if (err)
    {
        std::cerr << "hwloc: Failed to load topology" << std::endl;
        return 1;
    }
    err = hwloc_topology_export_xml(topology, filename, flags);
    if (err)
    {
        std::cerr << "hwloc: Failed to export xml" << std::endl;
        return 1;
    }

    return 0;
}

int execute(const std::string &command, const std::string &arg1, const std::string &arg2,
            const std::string &arg3, const std::string &arg4, const std::string &arg5, const std::string &arg6)
{
    pid_t pid = fork();
    if (pid == 0)
    {
        // Execute command through cli
        if (arg2 == "")
            execlp(command.c_str(), command.c_str(), arg1.c_str(), nullptr);
        else if (arg4 == "")
            execlp(command.c_str(), command.c_str(), arg1.c_str(), arg2.c_str(), arg3.c_str(), nullptr);
        else if (arg5 == "")
            execlp(command.c_str(), command.c_str(), arg1.c_str(), arg2.c_str(), arg3.c_str(), arg4.c_str(), nullptr);
        else
            execlp(command.c_str(), command.c_str(), arg1.c_str(), arg2.c_str(), arg3.c_str(), arg4.c_str(), arg5.c_str(), arg6.c_str(), nullptr);

        // Only executed when execlp failed
        std::cerr << "Failed executing: " << command << std::endl;
        exit(1);
    }
    else if (pid > 0)
    {
        int status;
        waitpid(pid, &status, 0);
        if (WIFEXITED(status) && WEXITSTATUS(status) != 0)
        {
            // Command execution failed
            std::cerr << "Failed to execute command: " << command << std::endl;
            return 1;
        }
    }
    else
    {
        std::cerr << "Fork Failed" << std::endl;
        return 1;
    }

    return 0;
}

void run_stream(const std::vector<int> &cpu, const int array_size, const int cache_line_size, const std::string &hostname,
                const std::string &sId, const int core, const double idle_power, const std::string &result)
{
    // Set affinity to current cpu inside this component
    change_affinity(cpu);

    set_omp(cpu);

    // Compile benchmark
    if (execute("make", "stream", "N=" + std::to_string(array_size), "D=" + hostname, "", "", "") != 0)
    {
        std::cerr << "Failed to compile benchmark file in " << std::string(hostname) << std::endl;
        return;
    }

    // Execute benchmark

    if (execute("./stream" + hostname, hostname, sId, std::to_string(core).c_str(), std::to_string(idle_power), std::to_string(cache_line_size), result) != 0)
    {
        std::cerr << "Failed to run benchmark in " << std::string(hostname) << ". Is the executable there?" << std::endl;
        return;
    }

    // Delete executable, so next iteration will have to create again.
    if (execute("rm", "stream" + hostname, "", "", "", "", "") != 0)
    {
        std::cerr << "Failed to delete the executable file stream" << std::string(hostname) << "." << std::endl;
        return;
    }
}

void run_dgemm(const std::vector<int> &cpu, const int array_size, const int limit, const std::string &hostname,
               const std::string &sId, const std::string &cId, const double idle_power, const std::string &result)
{
    // Set affinity to current cpu inside this component
    change_affinity(cpu);

    set_omp(cpu);

    // Compile benchmark
    if (execute("make", "dgemm", "N=" + std::to_string(array_size), "D=" + hostname, "", "", "") != 0)
    {
        std::cerr << "Failed to compile benchmark file in " << std::string(hostname) << std::endl;
        return;
    }

    // Execute benchmark
    if (execute("./dgemm" + hostname, hostname, sId, cId, std::to_string(idle_power), std::to_string(limit), result) != 0)
    {
        std::cerr << "Failed to run benchmark in " << std::string(hostname) << ". Is the executable there?" << std::endl;
        return;
    }

    // Delete executable, so next iteration will have to create again.
    if (execute("rm", "dgemm" + hostname, "", "", "", "", "") != 0)
    {
        std::cerr << "Failed to delete the executable file dgemm" << std::string(hostname) << "." << std::endl;
        return;
    }
}

/*
for memory:
thread = number of threads
limit = line size

for cpu:
thread = threadId,number of threads
limit = limit multiplier
*/
void run_nas(const std::vector<int> &cpu, const std::string &level, const std::string &node, const std::string &socket, const std::string &threads, const int &line_size, const std::string benchmark)
{
    std::string filename = "result/data/ranker_" + level + "_" + benchmark + "_" + node + ".csv";

    // Set affinity
    change_affinity(cpu);

    set_omp(cpu);

    std::string events = "-e LLC-load-misses,LLC-store-misses,fp_arith_inst_retired.scalar_double";

    std::string cmd;
    if (level == "core")
    { // Only FLOP, because there is no core level in memory
        events = "-e fp_arith_inst_retired.scalar_double";

        std::string cpuId = threads.substr(0, threads.find(","));
        for (char &c : cpuId)
        {
            if (c == ';')
                c = ',';
        }

        if (benchmark == "ep")
            cmd = "perf stat -a -C " + cpuId + " " + events + " ./NPB3.4.2/NPB3.4-OMP/bin/" + benchmark + ".A.x 2> " + node + "_tmp.txt";
        else
            cmd = "perf stat -a -C " + cpuId + " " + events + " ./NPB3.4.2/NPB3.4-OMP/bin/" + benchmark + ".A.x 2> " + node + "_tmp.txt";
    }
    else if (level == "socket")
    {
        if (benchmark == "ep")
            cmd = "perf stat --all-cpus --per-socket " + events + " ./NPB3.4.2/NPB3.4-OMP/bin/" + benchmark + ".C.x 2> " + node + "_tmp.txt";
        else
            cmd = "perf stat --all-cpus --per-socket " + events + " ./NPB3.4.2/NPB3.4-OMP/bin/" + benchmark + ".C.x 2> " + node + "_tmp.txt";
    }
    else
    {
        if (benchmark == "ep")
            cmd = "perf stat --all-cpus " + events + " ./NPB3.4.2/NPB3.4-OMP/bin/" + benchmark + ".D.x 2> " + node + "_tmp.txt";
        else
            cmd = "perf stat --all-cpus " + events + " ./NPB3.4.2/NPB3.4-OMP/bin/" + benchmark + ".D.x 2> " + node + "_tmp.txt";
    }

    if (execute("sh", "-c", cmd, "", "", "", "") != 0)
    {
        std::cerr << "Failed to run " + benchmark + " benchmark under perf in " << node << std::endl;
        return;
    }

    if (level == "core")
        write_nas_core_result(filename, node + "_tmp.txt", node, socket, threads);
    else
        write_nas_result(filename, node + "_tmp.txt", node, socket, threads, line_size);

    return;
}

void run(const std::string &benchmark, const std::string &level, const std::string &result, const int limit)
{
    uint32_t array_size = 0;

    char hostname[MPI_MAX_PROCESSOR_NAME];
    int len;
    MPI_Get_processor_name(hostname, &len);

    std::string filename = "result/topo/topology_" + std::string(hostname) + ".xml";

    init_topo(filename.c_str());

    // Root of sys-sage topology
    Topology *topo = new Topology();

    // Insert new node to the Topology
    Node *n = new Node(topo, 0);

    std::cout << "-- Parsing Hwloc output from " << filename << std::endl;

    // Parse hwloc to Topology
    if (parseHwlocOutput(n, filename) != 0)
    {
        std::cerr << "Failed to parse Hwloc output from " << filename << std::endl;
        return;
    }

    // Determine the benchmark level
    int ss_comp = SYS_SAGE_COMPONENT_NONE;
    if (level == "core")
        ss_comp = SYS_SAGE_COMPONENT_CORE;
    else if (level == "socket")
        ss_comp = SYS_SAGE_COMPONENT_CHIP;
    else
        ss_comp = SYS_SAGE_COMPONENT_NODE;

    // Get all subcomponents by the level
    std::vector<Component *> *comps = new std::vector<Component *>();
    topo->FindAllSubcomponentsByType(comps, ss_comp);

    // SocketID
    std::string sId = "-";

    int cache_line_size;

    std::string result_hostname = result + "_" + std::string(hostname);

    bool idle_power_recorded = false;
    std::vector<double> idle;

    // Iterate through all the subcomponents
    for (Component *comp : *comps)
    { // Setting array size and also skip iteration if chip is not socket.
        if (ss_comp == SYS_SAGE_COMPONENT_CORE)
        {
            Chip *chip = static_cast<Chip *>(comp->FindParentByType(SYS_SAGE_COMPONENT_CHIP));

            if (chip->GetName() != "socket")
                continue;

            sId = std::to_string(chip->GetId());

            std::vector<Component *> *cache = new std::vector<Component *>();
            chip->GetSubcomponentsByType(cache, SYS_SAGE_COMPONENT_CACHE);

            cache_line_size = static_cast<Cache *>(cache->front())->GetCacheLineSize();
            // For core benchmarking use half the capacity
            array_size = static_cast<Cache *>(cache->front())->GetCacheSize() / 8;
        }
        else if (ss_comp == SYS_SAGE_COMPONENT_CHIP)
        {
            if (comp->GetName() != "socket")
                continue;

            sId = std::to_string(comp->GetId());

            std::vector<Component *> *cache = new std::vector<Component *>();
            comp->GetSubcomponentsByType(cache, SYS_SAGE_COMPONENT_CACHE);
            cache_line_size = static_cast<Cache *>(cache->front())->GetCacheLineSize();

            // Number of double fit in cache
            array_size = static_cast<Cache *>(cache->front())->GetCacheSize() / 8;
        }
        else
        {
            std::vector<Component *> *chip = new std::vector<Component *>();
            comp->FindAllSubcomponentsByType(chip, SYS_SAGE_COMPONENT_CHIP);

            for (Component *c : *chip)
            {
                // If chip is socket, add L3 cache size to array_size, if not remove from comp so its threads doesn't get count
                if (c->GetName() == "socket")
                {
                    std::vector<Component *> *cache = new std::vector<Component *>();
                    c->GetSubcomponentsByType(cache, SYS_SAGE_COMPONENT_CACHE);

                    cache_line_size = static_cast<Cache *>(cache->front())->GetCacheLineSize();
                    array_size += static_cast<Cache *>(cache->front())->GetCacheSize() / 8;
                }
                else
                    comp->RemoveChild(c);
            }

            delete chip;
        }

        // Get the thread ids for setting affinity
        std::vector<Component *> *threads = new std::vector<Component *>();
        comp->FindAllSubcomponentsByType(threads, SYS_SAGE_COMPONENT_THREAD);

        // There is no cpu threads
        if (threads->empty())
            continue;

         if (benchmark == "stream")
        {
            // Get Idle power usage
            if (!idle_power_recorded)
            {
                record_power(result_hostname + "_idle", "w");
                sleep(2);
                record_power(result_hostname + "_idle", "a");

                // pkg0_watts, dram0_watts, pkg1_watts, dram1_watts, ...
                idle = parse_poll_power_idle(result_hostname + "_idle");

                idle_power_recorded = true;
            }
            // calculate idle power
            double idle_power = 0;
            if (sId == "-")
            {
                for (size_t i = 1; i < idle.size(); i = i + 2)
                    idle_power += idle[i];
            }
            else
            {
                idle_power = idle[stoi(sId) * 2 + 1];
            }

            for (size_t i = 2; i <= threads->size(); i = i + 2)
            {
                // CPUs to be tested.
                std::vector<int> cpu(i);

                std::transform(threads->begin(), threads->begin() + i, cpu.begin(),
                               [](Component *thread)
                               { return thread->GetId(); });

                run_stream(cpu, array_size * 4, cache_line_size, std::string(hostname), sId, i, idle_power, result_hostname + "_" + sId + ".csv");
            }

            // Get the Dynamic energy per line with least square
            calculate_dram_dynamic(result_hostname + "_" + sId + ".csv", result_hostname + ".csv");
            // CPUs to be tested.
            std::vector<int> cpu(threads->size());

            std::transform(threads->begin(), threads->end(), cpu.begin(),
                           [](Component *thread)
                           { return thread->GetId(); });

            run_nas(cpu, level, std::string(hostname), sId, "-," + to_string(threads->size()), cache_line_size, "ep");

            run_nas(cpu, level, std::string(hostname), sId, "-," + to_string(threads->size()), cache_line_size, "bt");

            run_nas(cpu, level, std::string(hostname), sId, "-," + to_string(threads->size()), cache_line_size, "sp");
        }
        else
        {
            //  Get tID
            std::string tId = "";

            if (level == "core")
            {
                for (size_t i = 0; i < threads->size(); i++)
                {
                    tId.append(std::to_string((*threads)[i]->GetId()));

                    if (i != threads->size() - 1)
                        tId.append(";");
                }
            }
            else
                tId = "-";

            // Get the power multiplier
            int start = result.find("-");
            std::string end = result.substr(start + 1, result.size());

            end.replace(1, 1, ".");

            double end_double = stod(end);
            // Core level start with 20% of the socket power limit
            if (level == "core")
            {
                variorum_cap_each_socket_power_limit(limit * end_double * 0.2);
                std::cout << "Socket power limit set to " << static_cast<int>(limit * end_double * 0.2) << "W" << std::endl;
            }
            else
            {
                variorum_cap_each_socket_power_limit(limit * end_double);
                std::cout << "Socket power limit set to " << static_cast<int>(limit * end_double) << "W" << std::endl;
            }

            // Get Idle power usage
            if (!idle_power_recorded)
            {

                record_power(result_hostname + "_idle", "w");
                sleep(2);
                record_power(result_hostname + "_idle", "a");

                // pkg0_watts, dram0_watts, pkg1_watts, dram1_watts, ...
                idle = parse_poll_power_idle(result_hostname + "_idle");

                idle_power_recorded = true;
            }

            double idle_power = 0;

            array_size = static_cast<int>(sqrt(array_size));

            if (sId == "-")
            {
                for (size_t i = 0; i < idle.size(); i = i + 2)
                    idle_power += idle[i];
            }
            else
            {
                idle_power = idle[stoi(sId) * 2];
            }

            // CPUs to be tested.
            std::vector<int> cpu(threads->size());

            std::transform(threads->begin(), threads->end(), cpu.begin(),
                           [](Component *thread)
                           { return thread->GetId(); });

             if (level == "core")
             {
                 run_dgemm(cpu, array_size, limit * end_double * 0.2, std::string(hostname), sId, tId, idle_power, result_hostname + ".csv");
             }
             else
             {
                 run_dgemm(cpu, array_size, limit * end_double, std::string(hostname), sId, tId, idle_power, result_hostname + ".csv");
             }

            if (end_double == 1.0 && level == "core")
            {
                run_nas(cpu, level, std::string(hostname), sId, tId + "," + to_string(threads->size()), 0, "ep");
            }

            variorum_cap_each_socket_power_limit(limit);
            std::cout << "Socket power limit set back to " << limit << "W" << std::endl;
        }
         delete threads;
    }

    delete comps;

    int mpi_size, mpi_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    // Get all hostnames
    char all_hostnames[mpi_size][MPI_MAX_PROCESSOR_NAME];
    MPI_Gather(hostname, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, all_hostnames, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, 0, MPI_COMM_WORLD);

    // Sort and write the data
    if (mpi_rank == 0)
    {
        std::vector<Row> final_data;
        std::vector<Row> ep_data;
        std::vector<Row> bt_data;
        std::vector<Row> sp_data;
        std::vector<Row> nolimit_data;

        for (char *c : all_hostnames)
        {
            to_vector(result + "_" + std::string(c) + ".csv", final_data);
            to_vector("result/data/ranker_" + level + "_ep_" + std::string(c) + ".csv", ep_data);
            to_vector("result/data/ranker_" + level + "_bt_" + std::string(c) + ".csv", bt_data);
            to_vector("result/data/ranker_" + level + "_sp_" + std::string(c) + ".csv", sp_data);

            // Get the peak performance for util
            if (benchmark == "dgemm")
                to_vector("result/data/ranker_" + level + "_dgemm-1,0_" + std::string(c) + ".csv", nolimit_data);
            else if (benchmark == "stream")
                to_vector(result + "_" + std::string(c) + ".csv", nolimit_data);
        }

        // Calculate and rank based on throughput

        std::string header_e;
        std::string header_t;

        if (benchmark == "dgemm")
        {
            header_e = "NodeId,SocketId,CPUId,NumberOfThreads,PowerLimit(W),FLOP,uJ/FLOP,Score\n";
            header_t = "NodeId,SocketId,CPUId,NumberOfThreads,PowerLimit(W),FLOPS(Job),PeakFLOPS,Utility,Throughput(FLOPS),Score\n";
        }
        else if (benchmark == "stream")
        {
            header_e = "NodeId,SocketId,LLCMisses,DynamicEnergyPerLine(uJ/line),Score\n";
            header_t = "NodeId,SocketId,PowerLimit(W),Bandwidth(Job),PeakBandwidth,Utility,Throughput(MB/s),Score\n";
        }

        for (int i = 0; i < 3; i++)
        {

            std::string suffix_e;
            std::string suffix_t;

            std::vector<std::vector<Row>> throughput_data;
            std::vector<Row> efficiency_data;

            if (i == 0)
            {
                suffix_e = "_ep_energy.csv";
                suffix_t = "_ep_throughput.csv";
                calculate_score(efficiency_data, throughput_data, final_data, nolimit_data, ep_data, benchmark);
            }
            else if (i == 1)
            {
                suffix_e = "_bt_energy.csv";
                suffix_t = "_bt_throughput.csv";
                calculate_score(efficiency_data, throughput_data, final_data, nolimit_data, bt_data, benchmark);
            }
            else
            {
                suffix_e = "_sp_energy.csv";
                suffix_t = "_sp_throughput.csv";
                calculate_score(efficiency_data, throughput_data, final_data, nolimit_data, sp_data, benchmark);
            }

            // Score based on througput
            for (size_t j = 0; j < throughput_data.size(); j++)
            {
                std::string path = result + (benchmark == "stream" ? to_string(j + 1) : "") + suffix_t;
                std::ofstream out_stream_t(path, ios::out | ios::trunc);
                if (!out_stream_t)
                {
                    std::cerr << "Error opening file" << std::endl;
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }

                out_stream_t << header_t;
                for (const Row r : throughput_data[j])
                {
                    for (size_t i = 0; i < r.columns.size(); i++)
                    {
                        out_stream_t << r.columns[i];

                        if (i < r.columns.size() - 1)
                            out_stream_t << ",";
                    }

                    out_stream_t << "\n";
                }

                out_stream_t.close();

                std::cout << benchmark << " result saved in " << path << std::endl;
            }
            // Score based on efficiency
            std::ofstream out_stream(result + suffix_e, ios::out | ios::trunc);
            if (!out_stream)
            {
                std::cerr << "Error opening file" << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
            }

            out_stream << header_e;
            for (const Row r : efficiency_data)
            {
                for (size_t i = 0; i < r.columns.size(); i++)
                {
                    out_stream << r.columns[i];

                    if (i < r.columns.size() - 1)
                        out_stream << ",";
                }

                out_stream << "\n";
            }

            out_stream.close();

            std::cout << benchmark << " result saved in " << result << suffix_e << std::endl;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
}
