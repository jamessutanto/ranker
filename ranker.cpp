#include "ranker.h"

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

extern "C"
{
#include "variorum.h"
}

#include "sys-sage.hpp"



void run(const std::string &, const std::string &, const std::string &);

void rank_dgemm_core(const int freqLimit, const std::string &result)
{
    MPI_Init(NULL, NULL);

    variorum_cap_each_core_frequency_limit(freqLimit);
    run("dgemm", "core", result);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}

void rank_dgemm_socket(const int powerLimit, const std::string &result)
{
    MPI_Init(NULL, NULL);

    variorum_cap_each_socket_power_limit(powerLimit);
    run("dgemm", "socket", result);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}

void rank_dgemm_node(const int powerLimit, const std::string &result)
{
    MPI_Init(NULL, NULL);

    variorum_cap_best_effort_node_power_limit(powerLimit);
    run("dgemm", "node", result);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}

void rank_stream_core(const int powerLimit, const std::string &result)
{
    MPI_Init(NULL, NULL);

    run("stream", "core", result);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}

void rank_stream_socket(const int powerLimit, const std::string &result)
{
    MPI_Init(NULL, NULL);

    run("stream", "socket", result);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}

void rank_stream_node(const int powerLimit, const std::string &result)
{
    MPI_Init(NULL, NULL);

    run("stream", "node", result);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
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

int execute(const std::string &command, const std::string &arg1, const std::string &arg2, const std::string &arg3)
{
    pid_t pid = fork();
    if (pid == 0)
    {
        // Execute command through cli
        if (arg2 == "")
            execlp(command.c_str(), command.c_str(), arg1.c_str(), nullptr);
        else if (arg3 == "")
            execlp(command.c_str(), command.c_str(), arg1.c_str(), arg2.c_str(), nullptr);
        else
            execlp(command.c_str(), command.c_str(), arg1.c_str(), arg2.c_str(), arg3.c_str(), nullptr);

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

void run(const std::string& benchmark, const std::string &level, const std::string &result)
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

    // Iterate through all the subcomponents
    for (Component *comp : *comps)
    { // Setting array size and also skip iteration if chip is not socket.
        if (ss_comp == SYS_SAGE_COMPONENT_CORE)
        {
            Chip *chip = static_cast<Chip *>(comp->FindParentByType(SYS_SAGE_COMPONENT_CHIP));

            if (chip->GetName() != "socket")
                continue;

            // L3 Cache size
            array_size = static_cast<Cache *>(chip->GetChildByType(SYS_SAGE_COMPONENT_CACHE))->GetCacheSize();
        }
        else if (ss_comp == SYS_SAGE_COMPONENT_CHIP)
        {
            if (comp->GetName() != "socket")
                continue;

            // L3 Cache size
            array_size = static_cast<Cache *>(comp->GetChildByType(SYS_SAGE_COMPONENT_CACHE))->GetCacheSize();
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
                    array_size += static_cast<Cache *>(socket->GetChildByType(SYS_SAGE_COMPONENT_CACHE))->GetCacheSize();
                }
                else
                    comp->RemoveChild(c);
            }

            delete chip;
        }

        // Get the thread ids for setting affinity
        std::vector<Component *> *threads = new std::vector<Component *>();
        comp->FindAllSubcomponentsByType(threads, SYS_SAGE_COMPONENT_THREAD);

        // CPUs to be tested.
        std::vector<int> cpu(threads->size());

        // There is no cpu threads
        if (threads->empty())
            continue;

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

        // std::cout << oss.str() << std::endl;

        // Set OMP to use certain threads only
        std::string ompPlaces = oss.str();
        setenv("OMP_PLACES", ompPlaces.c_str(), 1);

        // Set the number of threads for OMP to use
        std::string ompNumThreads = std::to_string(threads->size());
        setenv("OMP_NUM_THREADS", ompNumThreads.c_str(), 1);

        if (benchmark == "dgemm") {
            array_size = array_size/1000000;

        } else if (benchmark == "stream") {
            array_size = array_size*4;
        }

        // Compile benchmark
            if (execute("make", benchmark, "N=" + std::to_string(array_size), "D=" + std::string(hostname)) != 0)
            {
                std::cerr << "Failed to compile benchmark file." << std::endl;
                return;
            }

            // Execute benchmark
            if (execute("./" + benchmark + std::string(hostname), level, std::to_string(comp->GetId()), result) != 0)
            {
                std::cerr << "Failed to run benchmark. Is the executable there?" << std::endl;
                return;
            }

            // Delete executable, so next iteration will have to create again.
            if (execute("rm", benchmark + std::string(hostname), "", "") != 0)
            {
                std::cerr << "Failed to delete the executable file " << std::string(hostname) << "." << std::endl;
                return;
            }


        delete threads;
    }

    delete comps;

}
