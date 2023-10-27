#include <iostream>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <fcntl.h>
#include <cstdio>
#include <string>
#include <sstream>
#include "ranker.h"
#include <fstream>

std::string executeAndGetOutput() {
    char buffer[128];
    std::string result = "";
    
    FILE* pipe = popen("perf stat -e LLC-load-misses,LLC-store-misses ./NPB3.4.2/NPB3.4-OMP/bin/ep.D.x 2>&1", "r");
    if (!pipe) {
        return "popen failed!";
    }

    while (!feof(pipe)) {
        if (fgets(buffer, 128, pipe) != NULL)
            result += buffer;
    }

    pclose(pipe);
    return result;
}

long long getLLCLoadMisses(const std::string& perfOutput) {
    std::istringstream stream(perfOutput);
    std::string line;
    std::string numberStr;

    while (std::getline(stream, line)) {
        if (line.find("seconds time elapsed") != std::string::npos) {
            std::istringstream lineStream(line);

            long long misses;
            lineStream >> misses;
            return misses;
        }
    }
    return -1;  // Return -1 if not found
}

int main() {
    
/*     std::string output = executeAndGetOutput();

    long long llcMisses = getLLCLoadMisses(output);
    if (llcMisses != -1) {
        std::cout << "LLC-load-misses: " << llcMisses << std::endl;
    } else {
        std::cout << "Failed to find LLC-load-misses in perf output." << std::endl;
    }
    
    std::cout << output << std::endl; */
    rank_node();
    /* std::string perfOutput = "";

    FILE *pipe = popen("perf stat -e LLC-load-misses ./stream 2> testtest.txt", "r");
    if (!pipe)
    {
        std::cerr << "popen failed! fail to execute EP under perf" << std::endl;
        return -1;
    }
    std::cout << "after popen" << std::endl;

    char buffer[128];
    while (fgets(buffer, 128, pipe) != NULL)
    {
        std::cout << buffer;
        perfOutput += buffer;
    }

    pclose(pipe); */

/*     std::string ep_filename;

    std::string header;
    std::string events;

    std::string component ="memory";
    std::string level = "node";
    std::string node ="dsa";
    std::string limit = "2";

    if (component == "memory")
    {
        ep_filename = "ranker_" + level + "_ep_" + component + "_" + node + ".csv";

        header = "NodeId,SocketId,NumberOfThreads,LineSize(B),Time(S),LLCMisses,Bandwidth(MB/s)\n";
        events = "-e LLC-load-misses,LLC-store-misses";
    }
    else if (component == "cpu")
    {
        ep_filename = "ranker_" + level + "_ep_" + component + "-" + limit + "_" + node + ".csv";

        header = "NodeId,SocketId,ThreadId,NumberOfThreads,Time(S),TotalFLOP,FLOPS\n";
        events = "-e fp_arith_inst_retired.scalar_double";
    }

    std::string cmd;
    if (level == "core")
    {
        cmd = "perf stat " + events + " ./NPB3.4.2/NPB3.4-OMP/bin/ep.A.x > 2>&1";
    }
    else
    {
        cmd = "perf stat " + events + " ./NPB3.4.2/NPB3.4-OMP/bin/ep.A.x 2>&1";
    }
    std::cout << "before popen" << std::endl;
    std::cout << cmd << std::endl;

    std::string perfOutput = "";

    FILE *pipe = popen(cmd.c_str(), "r");
    if (!pipe)
    {
        std::cerr << "popen failed! fail to execute EP under perf" << std::endl;
        return -1;
    }
    std::cout << "after popen" << std::endl;

    char buffer[128];
    while (fgets(buffer, 128, pipe) != NULL)
    {
        std::cout << buffer;
        perfOutput += buffer;
    }

    pclose(pipe); */
/*     std::string a,b,c,d,e;
    a="perf";
    b="stat";
    c="-e fp_arith_inst_retired.scalar_double";
    d="./NPB3.4.2/NPB3.4-OMP/bin/ep.A.x";
    e="2>&1";
    execlp(a.c_str(), a.c_str(), b.c_str(), c.c_str(), d.c_str(), e.c_str(), nullptr); */
    /* std::ifstream infile("a_tmp.txt");
    std::string line;
    std::vector<double> llc_load_misses_values;

    while (std::getline(infile, line)) {
        if (line.find("LLC-load-misses") != std::string::npos) {
            double value;
            char junk;
            int ret = sscanf(line.c_str(), "S%*d %*d %lf %c", &value, &junk);
            
            // Check if sscanf successfully read the value and the next char is not a digit (it should be a space or a dot)
            if (ret == 2 && !isdigit(junk)) {
                llc_load_misses_values.push_back(value);
            }
        }
    }

    for (double val : llc_load_misses_values) {
        std::cout << val << std::endl;
    } */

    return 0;
}