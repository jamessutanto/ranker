#include "ranker.h"
#include "variorum_parser.hpp"

#include <iostream>
#include <cstring>
#include <hwloc.h>
#include <string>
#include <vector>
#include <sstream>
#include <iterator>
#include <sched.h>
#include <sys/wait.h>
#include <mpi.h>
#include <fstream>
#include <filesystem>

extern "C"
{
#include "variorum.h"
}

#include "sys-sage.hpp"

void run(const std::string &, const std::string &, const std::string &, const int);

int get_socket_power_limit();

int get_node_power_limit();

void rank_dgemm_core(const std::string &result)
{
    MPI_Init(NULL, NULL);

    int powerLimit;

    int limit = get_socket_power_limit();
    std::cout << "Socket power limit is " << limit << std::endl;

    powerLimit = limit;
    run("dgemm", "core", result + "-1.0", powerLimit);

    powerLimit = static_cast<int>(limit * 0.7);
    variorum_cap_each_socket_power_limit(powerLimit);
    std::cout << "Socket power limit set to " << powerLimit << std::endl;

    run("dgemm", "core", result + "-0.7", powerLimit);

    powerLimit = static_cast<int>(limit * 0.5);
    variorum_cap_each_socket_power_limit(powerLimit);
    std::cout << "Socket power limit set to " << powerLimit << std::endl;

    run("dgemm", "core", result + "-0.5", powerLimit);

    powerLimit = static_cast<int>(limit * 0.3);
    variorum_cap_each_socket_power_limit(powerLimit);
    std::cout << "Socket power limit set to " << powerLimit << std::endl;

    run("dgemm", "core", result + "-0.3", powerLimit);

    variorum_cap_each_socket_power_limit(limit);
    std::cout << "Socket power limit set back to " << limit << std::endl;

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}

void rank_dgemm_socket(const std::string &result)
{
    MPI_Init(NULL, NULL);

    int powerLimit;

    int limit = get_socket_power_limit();
    std::cout << "Socket power limit is " << limit << std::endl;

    powerLimit = limit;
    run("dgemm", "socket", result + "-1.0", powerLimit);

    powerLimit = static_cast<int>(limit * 0.7);
    variorum_cap_each_socket_power_limit(powerLimit);
    std::cout << "Socket power limit set to " << powerLimit << std::endl;

    run("dgemm", "socket", result + "-0.7", powerLimit);

    powerLimit = static_cast<int>(limit * 0.5);
    variorum_cap_each_socket_power_limit(powerLimit);
    std::cout << "Socket power limit set to " << powerLimit << std::endl;

    run("dgemm", "socket", result + "-0.5", powerLimit);

    powerLimit = static_cast<int>(limit * 0.3);
    variorum_cap_each_socket_power_limit(powerLimit);
    std::cout << "Socket power limit set to " << powerLimit << std::endl;

    run("dgemm", "socket", result + "-0.3", powerLimit);

    variorum_cap_each_socket_power_limit(limit);
    std::cout << "Socket power limit set back to " << limit << std::endl;

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}

void rank_dgemm_node(const std::string &result)
{
    MPI_Init(NULL, NULL);

    int powerLimit;

    int limit = get_node_power_limit();
    std::cout << "node power limit is " << limit << std::endl;

    powerLimit = limit;
    run("dgemm", "node", result + "-1", powerLimit);

    powerLimit = static_cast<int>(limit * 0.7);
    variorum_cap_best_effort_node_power_limit(powerLimit);
    std::cout << "node power limit set to " << powerLimit << std::endl;

    run("dgemm", "node", result + "-0.7", powerLimit);

    powerLimit = static_cast<int>(limit * 0.5);
    variorum_cap_best_effort_node_power_limit(powerLimit);
    std::cout << "node power limit set to " << powerLimit << std::endl;

    run("dgemm", "node", result + "-0.5", powerLimit);

    powerLimit = static_cast<int>(limit * 0.3);
    variorum_cap_best_effort_node_power_limit(powerLimit);
    std::cout << "node power limit set to " << powerLimit << std::endl;

    run("dgemm", "node", result + "-0.3", powerLimit);

    variorum_cap_best_effort_node_power_limit(limit);
    std::cout << "node power limit set back to " << limit << std::endl;

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}

void rank_stream_socket(const std::string &result)
{
    MPI_Init(NULL, NULL);

    run("stream", "socket", result, -1);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}

void rank_stream_node(const std::string &result)
{
    MPI_Init(NULL, NULL);

    run("stream", "node", result, -1);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}

void rank_srandom_socket(const std::string &result)
{
    MPI_Init(NULL, NULL);

    run("srandom", "socket", result, -1);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}

void rank_srandom_node(const std::string &result)
{
    MPI_Init(NULL, NULL);

    run("srandom", "node", result, -1);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}

void calculate_rank(const std::string &) {}

void rank_core()
{
    std::filesystem::path cwd = std::filesystem::current_path();

    rank_dgemm_core((cwd / "result/ranker_core_dgemm").string());

    calculate_rank("core");
}

void rank_socket()
{
    std::filesystem::path cwd = std::filesystem::current_path();

    rank_stream_socket(cwd / "result/ranker_socket_stream");
    // rank_srandom_socket(cwd / "result/ranker_socket_srandom");
    // rank_dgemm_socket(cwd / "result/ranker_socket_dgemm");

    calculate_rank("socket");
}

void rank_node()
{
    std::filesystem::path cwd = std::filesystem::current_path();

    rank_stream_node(cwd / "result/ranker_node_stream");
    rank_srandom_node(cwd / "result/ranker_node_srandom");
    rank_dgemm_node(cwd / "result/ranker_node_dgemm");

    calculate_rank("node");
}

int get_socket_power_limit()
{
    char hostname[MPI_MAX_PROCESSOR_NAME];
    int len;
    MPI_Get_processor_name(hostname, &len);

    record_power(hostname, "w");

    return socket_power_limit(hostname);
}

int get_node_power_limit()
{
    char hostname[MPI_MAX_PROCESSOR_NAME];
    int len;
    MPI_Get_processor_name(hostname, &len);

    record_power(hostname, "w");

    return node_power_limit(hostname);
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

void change_affinity(const std::vector<int> &cpu)
{
    cpu_set_t cpuset;

    // Set affinity only to CPU newAff
    CPU_ZERO(&cpuset);

    for (int i : cpu)
    {
        CPU_SET(i, &cpuset);
    }

    if (sched_setaffinity(0, sizeof(cpuset), &cpuset) == -1)
    {
        std::cerr << "Failed to set CPU affinity" << std::endl;
        return;
    }
}

int execute(const std::string &command, const std::string &arg1, const std::string &arg2,
            const std::string &arg3, const std::string &arg4, const std::string &arg5)
{
    pid_t pid = fork();
    if (pid == 0)
    {
        // Execute command through cli
        if (arg2 == "")
            execlp(command.c_str(), command.c_str(), arg1.c_str(), nullptr);
        else if (arg4 == "")
            execlp(command.c_str(), command.c_str(), arg1.c_str(), arg2.c_str(), arg3.c_str(), nullptr);
        else
            execlp(command.c_str(), command.c_str(), arg1.c_str(), arg2.c_str(), arg3.c_str(), arg4.c_str(), arg5.c_str(), nullptr);

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

struct Row
{
    std::vector<std::string> columns;
};

std::istream &operator>>(std::istream &str, Row &data)
{
    data.columns.clear();

    std::string line;
    std::getline(str, line);
    std::stringstream stream(line);
    std::string cell;

    while (std::getline(stream, cell, ','))
        data.columns.push_back(cell);

    return str;
}

bool compare_rows(const Row &r1, const Row &r2)
{
    return std::stod(r1.columns.back()) < std::stod(r2.columns.back());
}

bool compare_rows_reverse(const Row &r1, const Row &r2)
{
    return std::stod(r1.columns.back()) > std::stod(r2.columns.back());
}

void to_vector(const std::string &filename, std::vector<Row> &local)
{
    std::ifstream file;
    file.open(filename, std::ios::in);
    if (!file.is_open())
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
        return;
    }

    Row header;
    file >> header;

    Row row;
    while (file >> row)
    {
        local.push_back(row);
    }

    file.close();
}

void set_omp(const std::vector<int> &cpu)
{
    // Transform vector to string for OMP_PLACES
    std::ostringstream oss;
    oss << "{";
    std::move(cpu.begin(), cpu.end() - 1, std::ostream_iterator<int>(oss, ", "));
    oss << cpu.back();
    oss << "}";

    // Set OMP to use certain threads only
    std::string ompPlaces = oss.str();
    setenv("OMP_PLACES", ompPlaces.c_str(), 1);

    // Set the number of threads for OMP to use
    std::string ompNumThreads = std::to_string(cpu.size());
    setenv("OMP_NUM_THREADS", ompNumThreads.c_str(), 1);
}

void run_stream(const std::vector<int> &cpu, const int array_size, const std::string &hostname,
                const std::string &sId, const int core, const double idle_power, const std::string &result)
{
    // Set affinity to current cpu inside this component
    change_affinity(cpu);

    set_omp(cpu);

    // Compile benchmark
    if (execute("make", "stream", "N=" + std::to_string(array_size), "D=" + hostname, "", "") != 0)
    {
        std::cerr << "Failed to compile benchmark file in " << std::string(hostname) << std::endl;
        return;
    }

    // Execute benchmark
    if (execute("./stream" + hostname, hostname, sId, std::to_string(core).c_str(), std::to_string(idle_power), result) != 0)
    {
        std::cerr << "Failed to run benchmark in " << std::string(hostname) << ". Is the executable there?" << std::endl;
        return;
    }

    // Delete executable, so next iteration will have to create again.
    if (execute("rm", "stream" + hostname, "", "", "", "") != 0)
    {
        std::cerr << "Failed to delete the executable file stream" << std::string(hostname) << "." << std::endl;
        return;
    }
}

void run_srandom(const std::vector<int> &cpu, const int array_size, const std::string &hostname,
                 const std::string &sId, const int core, const double idle_power, const std::string &result)
{
    // Set affinity to current cpu inside this component
    change_affinity(cpu);

    set_omp(cpu);

    // Compile benchmark
    if (execute("make", "srandom", "N=" + std::to_string(array_size), "D=" + hostname, "", "") != 0)
    {
        std::cerr << "Failed to compile benchmark file in " << std::string(hostname) << std::endl;
        return;
    }

    // Execute benchmark
    if (execute("./srandom" + hostname, hostname, sId, std::to_string(core).c_str(), std::to_string(idle_power), result) != 0)
    {
        std::cerr << "Failed to run benchmark in " << std::string(hostname) << ". Is the executable there?" << std::endl;
        return;
    }

    // Delete executable, so next iteration will have to create again.
    if (execute("rm", "srandom" + hostname, "", "", "", "") != 0)
    {
        std::cerr << "Failed to delete the executable file srandom" << std::string(hostname) << "." << std::endl;
        return;
    }
}

void run_dgemm(const std::vector<int> &cpu, const int array_size, const std::string &hostname,
               const std::string &sId, const std::string &cId, const double idle_power, const std::string &result)
{
    // Set affinity to current cpu inside this component
    change_affinity(cpu);

    set_omp(cpu);

    // Compile benchmark
    if (execute("make", "dgemm", "N=" + std::to_string(array_size), "D=" + hostname, "", "") != 0)
    {
        std::cerr << "Failed to compile benchmark file in " << std::string(hostname) << std::endl;
        return;
    }

    // Execute benchmark
    if (execute("./dgemm" + hostname, hostname, sId, cId, std::to_string(idle_power), result) != 0)
    {
        std::cerr << "Failed to run benchmark in " << std::string(hostname) << ". Is the executable there?" << std::endl;
        return;
    }

    // Delete executable, so next iteration will have to create again.
    if (execute("rm", "dgemm" + hostname, "", "", "", "") != 0)
    {
        std::cerr << "Failed to delete the executable file dgemm" << std::string(hostname) << "." << std::endl;
        return;
    }
}

double calculate_efficiency(const std::string &path, const std::string &final)
{
    return 0;
}

void run(const std::string &benchmark, const std::string &level, const std::string &result, const int limit)
{
    uint32_t array_size = 0;

    char hostname[MPI_MAX_PROCESSOR_NAME];
    int len;
    MPI_Get_processor_name(hostname, &len);

    std::string filename = "result/topology_" + std::string(hostname) + ".xml";

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

    // SocketID and ThreadId
    std::string sId = "-";
    std::string cId = "-";

    std::string result_hostname = result + "_" + std::string(hostname);

    // Get Idle power usage

    record_power(result_hostname.c_str(), "w");
    sleep(2);
    record_power(result_hostname.c_str(), "a");

    // pkg0_joules, dram0_watts, pkg1_joules, dram1_watts, ...
    std::vector<double> idle(parse_poll_power_idle(result_hostname));

    // Iterate through all the subcomponents
    for (Component *comp : *comps)
    { // Setting array size and also skip iteration if chip is not socket.
        if (ss_comp == SYS_SAGE_COMPONENT_CORE)
        {
            Chip *chip = static_cast<Chip *>(comp->FindParentByType(SYS_SAGE_COMPONENT_CHIP));

            if (chip->GetName() != "socket")
                continue;

            sId = std::to_string(chip->GetId());
            // L3 Cache size
            array_size = static_cast<Cache *>(chip->GetChildByType(SYS_SAGE_COMPONENT_CACHE))->GetCacheSize();
        }
        else if (ss_comp == SYS_SAGE_COMPONENT_CHIP)
        {
            if (comp->GetName() != "socket")
                continue;

            sId = std::to_string(comp->GetId());
            // L3 Cache size (Cache size in bytes *2/8 (size of double))
            array_size = static_cast<Cache *>(comp->GetChildByType(SYS_SAGE_COMPONENT_CACHE))->GetCacheSize() / 4;
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
                    Chip *socket = static_cast<Chip *>(c);
                    array_size += static_cast<Cache *>(socket->GetChildByType(SYS_SAGE_COMPONENT_CACHE))->GetCacheSize() / 4;
                }
                else
                    comp->RemoveChild(c);
            }

            delete chip;
        }

        // Get the thread ids for setting affinity
        std::vector<Component *> *threads = new std::vector<Component *>();
        comp->FindAllSubcomponentsByType(threads, SYS_SAGE_COMPONENT_THREAD);

        if (level == "core")
        {
            for (size_t i = 0; i < threads->size(); i++)
            {
                cId = std::to_string((*threads)[i]->GetId());

                if (i != threads->size() - 1)
                    cId.append(";");
            }
        }

        // CPUs to be tested.
        std::vector<int> cpu(threads->size());

        // There is no cpu threads
        if (threads->empty())
            continue;

        if (benchmark != "dgemm")
        {
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

            for (size_t i = 1; i <= threads->size(); i++)
            {
                // CPUs to be tested.
                std::vector<int> cpu(i);

                std::transform(threads->begin(), threads->begin() + i, cpu.begin(),
                               [](Component *thread)
                               { return thread->GetId(); });

                if (benchmark == "stream")
                    run_stream(cpu, array_size * 4, std::string(hostname), sId, i, idle_power, result_hostname + "_" + sId);
                else if (benchmark == "srandom")
                    run_srandom(cpu, array_size * 4, std::string(hostname), sId, i, idle_power, result_hostname + "_" + sId);
            }

            calculate_efficiency(result_hostname + "_" + sId, result_hostname);
        }
        else
        {
            double idle_power = 0;
            if (sId == "-")
            {
                for (size_t i = 0; i < idle.size(); i = i + 2)
                    idle_power += idle[i];

                array_size = 2048;
            }
            else
            {
                idle_power = idle[stoi(sId) * 2];
                if(cId == "-")
                    array_size = 1024;
                else
                    array_size = 512;
            }

            // CPUs to be tested.
            std::vector<int> cpu(threads->size());

            std::transform(threads->begin(), threads->end(), cpu.begin(),
                           [](Component *thread)
                           { return thread->GetId(); });

            run_dgemm(cpu, array_size, std::string(hostname), sId, cId, idle_power, result_hostname);
        }
        /*
                // Get the thread Id
                std::transform(threads->begin(), threads->end(), cpu.begin(),
                               [](Component *thread)
                               { return thread->GetId(); });

                // Set affinity to current cpu inside this component
                change_affinity(cpu);

                // Transform vector to string for OMP_PLACES
                std::ostringstream oss;
                oss << "{";
                std::move(cpu.begin(), cpu.end() - 1, std::ostream_iterator<int>(oss, ", "));
                oss << cpu.back();
                oss << "}";

                // Set OMP to use certain threads only
                std::string ompPlaces = oss.str();
                setenv("OMP_PLACES", ompPlaces.c_str(), 1);

                // Set the number of threads for OMP to use
                std::string ompNumThreads = std::to_string(threads->size());
                setenv("OMP_NUM_THREADS", ompNumThreads.c_str(), 1);

                if (benchmark == "dgemm")
                {
                    array_size = 1024;
                }
                else if (benchmark == "stream" || benchmark == "srandom")
                {
                    array_size = array_size * 4;
                    std::cout << array_size << std::endl;
                }

                // Compile benchmark
                if (execute("make", benchmark, "N=" + std::to_string(array_size), "D=" + std::string(hostname), "") != 0)
                {
                    std::cerr << "Failed to compile benchmark file." << std::endl;
                    return;
                }

                // Execute benchmark
                if (execute("./" + benchmark + std::string(hostname), std::string(hostname), sId, cId, result + "_" + std::string(hostname)) != 0)
                {
                    std::cerr << "Failed to run benchmark. Is the executable there?" << std::endl;
                    return;
                }

                // Delete executable, so next iteration will have to create again.
                if (execute("rm", benchmark + std::string(hostname), "", "", "") != 0)
                {
                    std::cerr << "Failed to delete the executable file " << std::string(hostname) << "." << std::endl;
                    return;
                } */

        delete threads;
    }

    delete comps;

    int mpi_size, mpi_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    char all_hostnames[mpi_size][MPI_MAX_PROCESSOR_NAME];
    MPI_Gather(hostname, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, all_hostnames, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, 0, MPI_COMM_WORLD);

    std::vector<Row> final_data;

    for (char *c : all_hostnames)
    {
        to_vector(result + "_" + level + "_" + benchmark + "_" + std::string(c), final_data);
    }

    // Sort and write the data
    if (mpi_rank == 0)
    {
        std::string header;

        if (benchmark == "dgemm")
        {
            std::sort(final_data.begin(), final_data.end(), compare_rows);
            header = "NodeId, SocketId, ThreadId, PowerLimit(W), Time(S), PkgEnergy(J)\n";
        }
        else if (benchmark == "stream" || benchmark == "srandom")
        {
            std::sort(final_data.begin(), final_data.end(), compare_rows_reverse);
            header = "NodeId, SocketId, ThreadId, MemoryUsage(MiB), Time(S), Bandwidth(MB/s), DramPowerUsage(W), AccessFrequency, DramEfficiency(MB/J)\n";
        }

        std::ofstream out_stream(result + "_" + level + "_" + benchmark, ios::out | ios::trunc);
        if (!out_stream)
        {
            std::cerr << "Error opening file" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        out_stream << header;
        for (const Row r : final_data)
        {
            for (size_t i = 0; i < 3; i++)
            {
                out_stream << r.columns[i] << ",";
            }

            if (limit != -1)
                out_stream << limit << ",";

            for (size_t i = 3; i < r.columns.size(); i++)
            {
                out_stream << r.columns[i];

                if (i < r.columns.size() - 1)
                    out_stream << ",";
            }

            out_stream << "\n";
        }

        out_stream.close();

        std::cout << benchmark << "result saved in " << result << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);
}
