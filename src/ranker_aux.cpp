#include "ranker_aux.hpp"

#include <vector>
#include <unistd.h>
#include <filesystem>
#include <mpi.h>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <sched.h>

extern "C"
{
#include "variorum.h"
}

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

bool compare_last_reversed(const Row &r1, const Row &r2)
{
    size_t i = r1.columns.size();
    return std::stod(r1.columns[i - 1]) > std::stod(r2.columns[i - 1]);
}

bool compare_secondlast(const Row &r1, const Row &r2)
{
    size_t i = r1.columns.size();
    return std::stod(r1.columns[i - 2]) < std::stod(r2.columns[i - 2]);
}

void to_vector_header(const std::string &filename, std::vector<Row> &local)
{
    std::ifstream file;
    file.open(filename, std::ios::in);
    if (!file.is_open())
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
        return;
    }

    Row row;
    while (file >> row)
    {
        local.push_back(row);
    }

    file.close();
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

double RMSE(std::vector<double> &observed, std::vector<double> &predicted)
{
    double sum = 0.0;
    for (size_t i = 0; i < observed.size(); ++i)
    {
        double diff = observed[i] - predicted[i];
        sum += diff * diff;
    }

    double mse = sum / observed.size();
    return sqrt(mse);
}

void write_nas_core_result(const std::string &path_to, const std::string &path_from, const std::string &node, const std::string &socket, const std::string &threads)
{
    std::ifstream inFile(path_from);

    if (!inFile)
    {
        std::cerr << "Failed to open file." << std::endl;
        return;
    }

    std::string line;

    long long flop = 0;
    double time = 0;

    while (std::getline(inFile, line))
    { // if we run on socket level, choose the corresponding socket only
        if (line.find("fp_arith_inst_retired.scalar_double") != std::string::npos)
        {

            std::string number_str;
            for (char c : line)
            {
                if (std::isdigit(c))
                {
                    number_str += c;
                }
            }
            flop = stoll(number_str);
        }
        else if (line.find("seconds time elapsed") != std::string::npos)
        {
            std::string number_str;

            for (char c : line)
            {
                if (std::isdigit(c))
                {
                    number_str += c;
                }
                else if (c == ',')
                {
                    number_str += ".";
                }
            }
            time = std::stod(number_str);
        }
    }
    inFile.close();

    std::string header = "NodeId,SocketId,CPUId,NumberOfThreads,Time(S),TotalFLOP,FLOPS\n";

    std::ifstream alreadyExists(path_to);
    if (alreadyExists.is_open())
        header = "";

    // Open the file in append mode
    std::ofstream outFile(path_to, std::ios::app);

    // Check if the file was opened successfully
    if (!outFile)
    {
        std::cerr << "Unable to open file for writing." << std::endl;
        return;
    }

    // Write to the end of the file
    outFile << header;
    outFile << node << ",";
    outFile << socket << ",";
    outFile << threads << ",";
    outFile << time << ",";
    outFile << flop << ",";
    outFile << static_cast<long long>(flop / time) << "\n";

    // Close the file
    outFile.close();
}

void write_nas_result(const std::string &path_to, const std::string &path_from, const std::string &node, const std::string &socket, const std::string &threads, const int &cache_line_size)
{
    std::ifstream inFile(path_from);

    if (!inFile)
    {
        std::cerr << "Failed to open file." << std::endl;
        return;
    }

    std::string line;

    long long access = 0;
    long long flop = 0;
    double time = 0;

    while (std::getline(inFile, line))
    { // if we run on socket level, choose the corresponding socket only
        if (line.find("LLC-store-misses") != std::string::npos)
        {
            if (socket == "-")
            {
                std::string number_str;
                for (char c : line)
                {
                    if (std::isdigit(c))
                    {
                        number_str += c;
                    }
                }
                access += stoll(number_str);
            }
            else if (line.find("S" + socket) != std::string::npos)
            { // Remove the dot.
                line.erase(std::remove(line.begin(), line.end(), '.'), line.end());

                long long value;
                char junk;
                int ret = sscanf(line.c_str(), "S%*d %*d %lld %c", &value, &junk);

                // Check if sscanf successfully read the value and the next char is not a digit
                if (ret == 2 && !isdigit(junk))
                {
                    access += value;
                }
            }
        }
        else if (line.find("LLC-load-misses") != std::string::npos)
        {
            if (socket == "-")
            {
                std::string number_str;
                for (char c : line)
                {
                    if (std::isdigit(c))
                    {
                        number_str += c;
                    }
                }
                access += stoll(number_str);
            }
            else if (line.find("S" + socket) != std::string::npos)
            { // Remove the dot.
                line.erase(std::remove(line.begin(), line.end(), '.'), line.end());

                long long value;
                char junk;
                int ret = sscanf(line.c_str(), "S%*d %*d %lld %c", &value, &junk);

                // Check if sscanf successfully read the value and the next char is not a digit
                if (ret == 2 && !isdigit(junk))
                {
                    access += value;
                }
            }
        }
        else if (line.find("fp_arith_inst_retired.scalar_double") != std::string::npos)
        { // Node
            if (socket == "-")
            {
                std::string number_str;
                for (char c : line)
                {
                    if (std::isdigit(c))
                    {
                        number_str += c;
                    }
                }
                flop = stoll(number_str);
            } // Socket
            else if (line.find("S" + socket) != std::string::npos)
            { // Remove the dot.
                line.erase(std::remove(line.begin(), line.end(), '.'), line.end());

                long long value;
                char junk;
                int ret = sscanf(line.c_str(), "S%*d %*d %lld %c", &value, &junk);

                // Check if sscanf successfully read the value and the next char is not a digit
                if (ret == 2 && !isdigit(junk))
                {
                    flop = value;
                }
            }
        }
        else if (line.find("seconds time elapsed") != std::string::npos)
        {
            std::string number_str;

            for (char c : line)
            {
                if (std::isdigit(c))
                {
                    number_str += c;
                }
                else if (c == ',')
                {
                    number_str += ".";
                }
            }
            time = std::stod(number_str);
        }
    }
    inFile.close();

    std::string header = "NodeId,SocketId,CPUId,NumberOfThreads,CacheLineSize(B),Time(S),LLCMisses,Bandwidth(MB/s),TotalFLOP,FLOPS\n";

    std::ifstream alreadyExists(path_to);
    if (alreadyExists.is_open())
        header = "";

    // Open the file in append mode
    std::ofstream outFile(path_to, std::ios::app);

    // Check if the file was opened successfully
    if (!outFile)
    {
        std::cerr << "Unable to open file for writing." << std::endl;
        return;
    }

    // Write to the end of the file
    outFile << header;
    outFile << node << ",";
    outFile << socket << ",";
    outFile << threads << ",";
    outFile << cache_line_size << ",";
    outFile << time << ",";
    outFile << access << ",";
    outFile << static_cast<double>(access * cache_line_size / 1000000) / time << ",";
    outFile << flop << ",";
    outFile << static_cast<u_long>(flop / time) << "\n";

    // Close the file
    outFile.close();
}

void write_stream_result(const std::string &path_to, const std::string &path_from, const std::string &node, const std::string &socket, const std::string &threads, const int &cache_line_size)
{
    std::ifstream inFile(path_from);

    if (!inFile)
    {
        std::cerr << "Failed to open file." << std::endl;
        return;
    }

    std::string line;

    long long access = 0;
    double energy = 0;
    double time = 0;

    while (std::getline(inFile, line))
    { // if we run on socket level, choose the corresponding socket only
        if (line.find("LLC-store-misses") != std::string::npos)
        {
            if (socket == "-")
            {
                std::string number_str;
                for (char c : line)
                {
                    if (std::isdigit(c))
                    {
                        number_str += c;
                    }
                }
                access += stoll(number_str);
            }
            else if (line.find("S" + socket) != std::string::npos)
            { // Remove the dot.
                line.erase(std::remove(line.begin(), line.end(), '.'), line.end());

                long long value;
                char junk;
                int ret = sscanf(line.c_str(), "S%*d %*d %lld %c", &value, &junk);

                // Check if sscanf successfully read the value and the next char is not a digit
                if (ret == 2 && !isdigit(junk))
                {
                    access += value;
                }
            }
        }
        else if (line.find("LLC-load-misses") != std::string::npos)
        {
            if (socket == "-")
            {
                std::string number_str;
                for (char c : line)
                {
                    if (std::isdigit(c))
                    {
                        number_str += c;
                    }
                }
                access += stoll(number_str);
            }
            else if (line.find("S" + socket) != std::string::npos)
            { // Remove the dot.
                line.erase(std::remove(line.begin(), line.end(), '.'), line.end());

                long long value;
                char junk;
                int ret = sscanf(line.c_str(), "S%*d %*d %lld %c", &value, &junk);

                // Check if sscanf successfully read the value and the next char is not a digit
                if (ret == 2 && !isdigit(junk))
                {
                    access += value;
                }
            }
        }
        else if (line.find("power/energy-ram/") != std::string::npos)
        {
            if (socket == "-")
            {
                std::string number_str;
                for (char c : line)
                {
                    if (std::isdigit(c))
                    {
                        number_str += c;
                    }
                    else if (c == ',')
                    {
                        number_str += ".";
                    }
                }
                energy = stod(number_str);
            }
            else if (line.find("S" + socket) != std::string::npos)
            { // Replace comma with dot to make it double
                for (char &c : line)
                {
                    if (c == ',')
                    {
                        c = '.';
                    }
                }

                double value;
                char junk;
                int ret = sscanf(line.c_str(), "S%*d %*d %lf %c", &value, &junk);

                // Check if sscanf successfully read the value and the next char is not a digit
                if (ret == 2 && !isdigit(junk))
                {
                    energy = value;
                }
            }
        }
        else if (line.find("seconds time elapsed") != std::string::npos)
        {
            std::string number_str;

            for (char c : line)
            {
                if (std::isdigit(c))
                {
                    number_str += c;
                }
                else if (c == ',')
                {
                    number_str += ".";
                }
            }
            time = std::stod(number_str);
        }
    }

    inFile.close();

    std::string header = "NodeId,SocketId,NumberOfThreads,LLCMisses,Time(S),Bandwidth(MB/s),EnergyUsage(J),PowerUsage(W),AccessFrequency(line/s)\n";

    std::ifstream alreadyExists(path_to);
    if (alreadyExists.is_open())
        header = "";

    // Open the file in append mode
    std::ofstream outFile(path_to, std::ios::app);

    // Check if the file was opened successfully
    if (!outFile)
    {
        std::cerr << "Unable to open file for writing." << std::endl;
        return;
    }

    // Write to the end of the file
    outFile << header;
    outFile << node << ",";
    outFile << socket << ",";
    outFile << threads << ",";
    outFile << access << ",";
    outFile << time << ",";
    outFile << static_cast<double>(access * cache_line_size / 1000000) / time << ",";
    outFile << energy << ",";
    outFile << energy / time << ",";
    outFile << static_cast<u_long>(access / time) << "\n";

    // Close the file
    outFile.close();
}

void calculate_score(std::vector<Row> &result_e, std::vector<std::vector<Row>> &result_t, const std::vector<Row> &data, const std::vector<Row> &data_nolimit, const std::vector<Row> &nas_data, const std::string &benchmark)
{
    /*
        data = "NodeId,SocketId,CPUId,NumberOfThreads,PowerLimit(W),IdleCPUPower(W),AvgTime(S),SDTime,AvgEnergyUsage(J),SDEnergy,AvgPowerConsumption(W),SDPowerConsumption,TotalFLOP,AvgFLOPS,SDFLOPS,AvguJ/FLOP,SDuJ/FLOP\n";
        nas = "NodeId,SocketId,CPUId,NumberOfThreads,CacheLineSize(B),Time(S),LLCMisses,Bandwidth(MB/s),TotalFLOP,FLOPS"
        throughput = "NodeId,SocketId,CPUId,NumberOfThreads,PowerLimit(W),Utility,Throughput(FLOPS),Score\n";
        energy = "NodeId,SocketId,CPUId,NumberOfThreads,PowerLimit(W),FLOP,uJ/FLOP,Score\n"; */

    if (benchmark == "dgemm")
    {
        // Find PeakFLOPS
        double max_flops = 0;
        for (size_t i = 0; i < data.size(); i++)
        {
            double flops = stod(data[i].columns[13]);
            if (max_flops < flops)
                max_flops = flops;
        }

        result_t.push_back(std::vector<Row>());

        for (size_t i = 0; i < data.size(); i++)
        { // Throughput-oriented
            Row r;

            if (data[i].columns[1] != nas_data[i].columns[1])
            {
                std::cerr << "not matched" << std::endl;
                return;
            }

            r.columns.push_back(data[i].columns[0]); // NodeId
            r.columns.push_back(data[i].columns[1]); // SocketId
            r.columns.push_back(data[i].columns[2]); // CPUId
            r.columns.push_back(data[i].columns[3]); // NumberOfThreads
            r.columns.push_back(data[i].columns[4]); // PowerLimit

            double flops_nas = stod(nas_data[i].columns[9]);
            r.columns.push_back(std::to_string(flops_nas)); // FLOPS of job

            double peak_flops = stod(data_nolimit[i].columns[13]);
            r.columns.push_back(std::to_string(peak_flops)); // PeakFLOPS of this component with no limit

            double util = flops_nas / peak_flops;
            r.columns.push_back(std::to_string(util)); // Util

            double throughput = stod(data[i].columns[13]); // FLOPS
            r.columns.push_back(std::to_string(throughput));

            r.columns.push_back(std::to_string(util * throughput / max_flops)); // Score

            result_t[0].push_back(r);

            // Energy-oriented
            Row m;
            m.columns.push_back(data[i].columns[0]); // NodeId
            m.columns.push_back(data[i].columns[1]); // SocketId
            m.columns.push_back(data[i].columns[2]); // CPUId
            m.columns.push_back(data[i].columns[3]); // NumberOfThreads
            m.columns.push_back(data[i].columns[4]); // PowerLimit

            double uj_flop = stod(data[i].columns[15]);
            double flop = stod(nas_data[i].columns[8]);

            m.columns.push_back(std::to_string(flop)); // FLOP

            m.columns.push_back(std::to_string(uj_flop)); // uJ/FLOP

            m.columns.push_back(std::to_string(flop * uj_flop * -1)); // Score

            result_e.push_back(m);
        }
        /*      data = "NodeId,SocketId,PeakBandwidth(MB/s),EstimatedMemoryStaticPower(W),DynamicEnergyPerLine(uJ/line),HighestPower(W),BandwidthPerPower((MB/s)/W),Constant";
                nas          = "NodeId,SocketId,CPUId,NumberOfThreads,CacheLineSize(B),Time(S),LLCMisses,Bandwidth(MB/s),TotalFLOP,FLOPS"
                throughput = "NodeId,SocketId,PowerLimit,BW,PeakBw,Utility,Throughput(MB/s),Score\n"; */
    }
    else if (benchmark == "stream")
    {
        // Find max bandwidth of all components
        double max_bw = 0;
        int lowest = stoi(data[0].columns[3]);
        int highest = stoi(data[0].columns[5]);

        std::vector<std::pair<double, double>> bw_p;
        for (size_t i = 0; i < data.size(); i++)
        {
            double bw = stod(data[i].columns[2]);
            if (max_bw < bw)
                max_bw = bw;

            int static_power = stoi(data[i].columns[3]);
            if (lowest > static_power)
                lowest = static_power;

            int max_power = stoi(data[i].columns[5]);
            if (highest < max_power)
                highest = max_power;

            bw_p.push_back(std::make_pair(stod(data[i].columns[6]), stod(data[i].columns[7])));
        }

        int interval = (highest - lowest) / 10;

        std::vector<double> max_bw_p;
        for (int i = 0; i < interval; i++)
        {
            double power = lowest + 10 * (i + 1);
            double max_bww = bw_p[0].first * power + bw_p[0].second;
            if(max_bww > stod(data[0].columns[2]))
                break;
            
            bool stop=false;

            for (size_t k = 1; k < data.size(); k++)
            {
                double d = bw_p[k].first * power + bw_p[k].second;

                if(d > stod(data[k].columns[2]))
                    stop=true;
                
                if (max_bww < d)
                    max_bww = d;
            }
            if(stop)
                break;
            
            max_bw_p.push_back(max_bww);

            // Initialize
            result_t.push_back(std::vector<Row>());
        }

        // For without limit;
        result_t.push_back(std::vector<Row>());

        for (size_t i = 0; i < data.size(); i++)
        {
            Row r;

            if (data[i].columns[1] != nas_data[i].columns[1])
            {
                std::cerr << "not matched" << std::endl;
                return;
            }

            r.columns.push_back(data[i].columns[0]); // NodeId
            r.columns.push_back(data[i].columns[1]); // SocketId
            r.columns.push_back("-"); //No Limit
            
            double bandwidth_nas = stod(nas_data[i].columns[7]);
            r.columns.push_back(std::to_string(bandwidth_nas)); // Bandwidht(J)

            double peak_bw = stod(data_nolimit[i].columns[2]);
            r.columns.push_back(std::to_string(peak_bw)); // Peak BW without limit

            double util = bandwidth_nas / peak_bw; // How to get dynamic energy/ line with different power limit(?)
            r.columns.push_back(std::to_string(util));

            double throughput = stod(data[i].columns[2]); // BW
            r.columns.push_back(std::to_string(throughput));

            r.columns.push_back(std::to_string(util * throughput / max_bw)); // Score

            result_t[0].push_back(r);

            // Rank with power limit based on BW vs Power
            for (size_t j = 1; j < result_t.size(); j++)
            {   
                int power = lowest + (10 * j);
                double tp = power * bw_p[i].first + bw_p[i].second;
                /* if (tp > peak_bw)
                {
                    break;
                } */
                
                Row s;
                s.columns.push_back(data[i].columns[0]); // NodeId
                s.columns.push_back(data[i].columns[1]); // SocketId

                
                s.columns.push_back(std::to_string(power)); // Power Limit

                s.columns.push_back(std::to_string(bandwidth_nas)); // Bandwidht(Job)
                s.columns.push_back(std::to_string(peak_bw));       // Peak BW without limit

                s.columns.push_back(std::to_string(util));

                s.columns.push_back(std::to_string(tp)); // Throughput


                s.columns.push_back(std::to_string(util * tp / max_bw_p[j - 1]));


                result_t[j].push_back(s);
            }

            // Energy
            Row m;
            m.columns.push_back(data[i].columns[0]); // NodeId
            m.columns.push_back(data[i].columns[1]); // SocketId

            double requests = stod(nas_data[i].columns[6]); // LLC Misses
            m.columns.push_back(std::to_string(requests));

            double uj_line = stod(data[i].columns[4]); // Dynamic Energy per Line
            m.columns.push_back(std::to_string(uj_line));

            m.columns.push_back(std::to_string(requests * uj_line * -1)); // Score

            result_e.push_back(m);
        }
    }

    for (size_t i = 0; i < result_t.size(); i++)
        std::sort(result_t[i].begin(), result_t[i].end(), compare_last_reversed);

    std::sort(result_e.begin(), result_e.end(), compare_last_reversed);
}

void calculate_efficiency(std::vector<Row> &result, const std::vector<Row> &data, const std::vector<Row> &data_nolimit, const std::vector<Row> &nas_data, const std::string &benchmark)
{
    /*
        data = "NodeId,SocketId,CPUId,NumberOfThreads,PowerLimit(W),IdleCPUPower(W),AvgTime(S),SDTime,AvgEnergyUsage(J),SDEnergy,AvgPowerConsumption(W),SDPowerConsumption,TotalFLOP,AvgFLOPS,SDFLOPS,AvguJ/FLOP,SDuJ/FLOP\n";
                nas = NodeId,SocketId,CPUId,NumberOfThreads,Time(S),TotalFLOP,FLOPS
        energy = "NodeId,SocketId,CPUId,NumberOfThreads,PowerLimit(W),FLOP,uJ/FLOP,Score\n"; */
    if (benchmark == "dgemm")
    {
        for (size_t i = 0; i < data.size(); i++)
        {
            double uj_flop = stod(data[i].columns[15]);
            double flop = stod(nas_data[i].columns[5]);

            Row r;

            if (data[i].columns[1] != nas_data[i].columns[1])
            {
                std::cerr << "not matched" << std::endl;
                return;
            }

            r.columns.push_back(data[i].columns[0]); // NodeId
            r.columns.push_back(data[i].columns[1]); // SocketId
            r.columns.push_back(data[i].columns[2]); // CPUId
            r.columns.push_back(data[i].columns[3]); // NumberOfThreads
            r.columns.push_back(data[i].columns[4]); // PowerLimit

            r.columns.push_back(std::to_string(flop)); // FLOP

            r.columns.push_back(std::to_string(uj_flop)); // uJ/FLOP

            r.columns.push_back(std::to_string(-flop * uj_flop * -1)); // Score

            result.push_back(r);
        }
        /*data = "NodeId,SocketId,EstimatedMemoryStaticPower(W),PeakBandwidth(MB/s),DynamicEnergyPerLine(uJ/line),RMSE\n"
        nas          = NodeId,SocketId,NumberOfThreads,LineSize(B),Time(S),LLCMisses,Bandwidth(MB/s)
        energy = "NodeId,SocketId,LLCMisses,DynamicEnergyPerLine(uJ/line),Score\n"; */
    }
    else if (benchmark == "stream")
    {
        for (size_t i = 0; i < data.size(); i++)
        {
            Row r;

            if (data[i].columns[1] != nas_data[i].columns[1])
            {
                std::cerr << "not matched" << std::endl;
                return;
            }

            r.columns.push_back(data[i].columns[0]); // NodeId
            r.columns.push_back(data[i].columns[1]); // SocketId

            double requests = stod(nas_data[i].columns[5]); // LLC Misses
            r.columns.push_back(std::to_string(requests));

            double uj_line = stod(data[i].columns[4]); // Dynamic Energy per Line
            r.columns.push_back(std::to_string(uj_line));

            r.columns.push_back(std::to_string(requests * uj_line * -1)); // Score

            result.push_back(r);
        }
    }

    std::sort(result.begin(), result.end(), compare_last_reversed);
}

void calculate_dram_dynamic(const std::string &path, const std::string &final)
{
    std::string header = "NodeId,SocketId,PeakBandwidth(MB/s),EstimatedMemoryStaticPower(W),DynamicEnergyPerLine(uJ/line),HighestPower(W),BandwidthPerPower((MB/s)/W),Constant\n";

    std::ofstream out_stream(final, std::ios::out | std::ios::app);
    if (!out_stream)
    {
        std::cerr << "Error opening file" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    out_stream.seekp(0, std::ios_base::end);
    if (out_stream.tellp() == 0)
        out_stream << header;

    std::vector<Row> data;
    to_vector_header(path, data);

    //"NodeId, SocketId, TotalMemoryRequest(MiB), NumberOfThreads, AvgTime(S), SDTime, Bandwidth(MB/s), SDBandwidth, AvgEnergyUsage(J), SDEnergyUsage, AvgPowerUsage(W), SDPowerUsage, AccessFrequency(line/s), SDAccessFrequency\n"
    Row head = data.front();
    size_t i_avg_power = 10;
    size_t i_access_freq = 12;
    size_t i_node_id = 0;
    size_t i_socket_id = 1;
    size_t i_bandwidth = 6;

    double peak_bandwidth = 0;

    double max_power = 0;

    std::vector<double> access_freq;
    std::vector<double> avg_power;
    std::vector<double> bw;
    for (size_t i = 1; i < data.size(); i++)
    {
        access_freq.push_back(stod(data[i].columns[i_access_freq]));

        double power = stod(data[i].columns[i_avg_power]);
        avg_power.push_back(power);

        if (max_power < power)
            max_power = power;

        double bandwidth = stod(data[i].columns[i_bandwidth]);
        bw.push_back(bandwidth);
        if (peak_bandwidth < bandwidth)
            peak_bandwidth = bandwidth;
    }

    std::pair<double, double> pair = linear_regression(access_freq, avg_power);

    std::vector<double> predicted;
    for (size_t i = 0; i < access_freq.size(); i++)
    {
        predicted.push_back(access_freq[i] * pair.first + pair.second);
    }

    double rmse = RMSE(avg_power, predicted);

    std::pair<double, double> bw_w = linear_regression(avg_power, bw);

    std::vector<double> predicted1;
    for (size_t i = 0; i < access_freq.size(); i++)
    {
        predicted1.push_back(avg_power[i] * bw_w.first + bw_w.second);
    }

    double rmse1 = RMSE(bw, predicted1);

    out_stream << data[1].columns[i_node_id] << ",";
    out_stream << data[1].columns[i_socket_id] << ",";
    out_stream << peak_bandwidth << ",";
    out_stream << pair.second << ",";
    out_stream << pair.first * 1000000 << ",";
    //out_stream << rmse << ",";
    out_stream << max_power << ",";
    out_stream << bw_w.first << ",";
    out_stream << bw_w.second << "\n";

    out_stream.close();
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

int socket_power_limit(const char *node)
{
    FILE *power = fopen(node, "r");
    fseek(power, 0, SEEK_SET);

    // 25 for powmon and timestamp, 30 for recording per metric
    int buffer_size = 25 + 30 * 5 * 4;
    char header[buffer_size];

    fgets(header, sizeof(header), power);
    std::string a = std::string(header);

    if (a.find("time", 6) != std::string::npos)
    {
        char row1[buffer_size];
        fgets(row1, sizeof(row1), power);
        a = std::string(row1);
    }

    fclose(power);

    int start = 0;
    int end = -1;

    size_t pos = 0;
    int cnt = 0;

    // Removing _POWMON and time
    pos = a.find(" ", pos);
    pos++;
    pos = a.find(" ", pos);

    while ((pos = a.find(" ", pos)) != std::string::npos)
    {
        cnt++;
        if (cnt == 2)
        {
            start = pos;
        }
        else if (cnt == 3)
        {
            end = pos;
            break;
        }
        pos++;
    }
    return stoi(a.substr(start, end - start));
}

int node_power_limit(const char *node)
{
    FILE *power = fopen(node, "r");
    fseek(power, 0, SEEK_SET);

    // 25 for powmon and timestamp, 30 for recording per metric
    int buffer_size = 25 + 30 * 5 * 4;
    char header[buffer_size];
    char row1[buffer_size];

    fgets(header, sizeof(header), power);
    fgets(row1, sizeof(row1), power);
    std::string a = std::string(row1);

    fclose(power);

    std::vector<int> start;
    std::vector<int> end;

    size_t pos = 0;
    int cnt = 0;

    // Removing _POWMON and time
    pos = a.find(" ", pos);
    pos++;
    pos = a.find(" ", pos);

    while ((pos = a.find(" ", pos)) != std::string::npos)
    {
        cnt++;
        if (cnt == 2)
        {
            start.push_back(pos);
        }
        else if (cnt == 3)
        {
            end.push_back(pos);
        }
        else if (cnt == 5)
        {
            cnt = 0;
        }
        pos++;
    }

    int limit = 0;
    for (size_t i = 0; i < start.size(); i++)
    {
        limit += stoi(a.substr(start[i], end[i] - start[i]));
    }

    return limit;
}

std::pair<double, double> parse_poll_power_pkg(const std::string &node, const std::string &socket)
{

    FILE *power = fopen(node.c_str(), "r");
    fseek(power, 0, SEEK_SET);

    // 25 for powmon and timestamp, 30 for recording per metric
    int buffer_size = 25 + 30 * 5 * 4;
    char header[buffer_size];
    char initial[buffer_size];
    char finish[buffer_size];

    fgets(header, sizeof(header), power);

    fgets(initial, sizeof(initial), power);

    std::string hd = std::string(header);
    std::string a;
    std::string b;

    // Header not printed
    if (hd.find("time", 6) == std::string::npos)
    {
        a = hd;
        b = std::string(initial);
    }
    else
    {
        a = std::string(initial);
        fgets(finish, sizeof(finish), power);
        b = std::string(finish);
    }
    fclose(power);

    // Position of pkg_joules
    std::vector<int> start, end;
    // Position of time
    int start_time1, end_time1, start_time2, end_time2;

    // Usage
    double usage = 0;
    double time = 0;

    size_t pos1 = 0;
    size_t pos2 = 0;
    int cnt = 0;

    // Get time
    pos1 = a.find(" ", pos1);
    start_time1 = pos1;
    pos1++;
    pos1 = a.find(" ", pos1);
    end_time1 = pos1;

    pos2 = b.find(" ", pos2);
    start_time2 = pos2;
    pos2++;
    pos2 = b.find(" ", pos2);
    end_time2 = pos2;

    while ((pos2 = b.find(" ", pos2)) != std::string::npos)
    {
        cnt++;
        if (cnt == 1)
        {
            start.push_back(pos2);
        }
        else if (cnt == 2)
        {
            end.push_back(pos2);
        }
        else if (cnt == 5)
        {
            cnt = 0;
        }
        pos2++;
    }

    // Sum usage of all sockets for node level
    if (socket == "-")
    {
        for (size_t i = 0; i < start.size(); i++)
        {
            usage += stod(b.substr(start[i], end[i] - start[i]));
        }
    }
    else
    { // Get individual usage
        int numS = stoi(socket);

        usage = stod(b.substr(start[numS], end[numS] - start[numS]));
    }

    time = stod(b.substr(start_time2, end_time2 - start_time2));
    time = time - stod(a.substr(start_time1, end_time1 - start_time1));
    time = time / 1000; // in second

    std::filesystem::remove(node);

    return std::make_pair(usage, usage / time);
}

std::pair<double, double> parse_poll_power_dram(const std::string &node, const std::string &socket)
{

    FILE *power = fopen(node.c_str(), "r");
    fseek(power, 0, SEEK_SET);

    // 25 for powmon and timestamp, 30 for recording per metric
    int buffer_size = 25 + 30 * 5 * 4;
    char header[buffer_size];
    char initial[buffer_size];
    char finish[buffer_size];

    fgets(header, sizeof(header), power);

    fgets(initial, sizeof(initial), power);

    std::string hd = std::string(header);
    std::string a;
    std::string b;

    // Header not printed
    if (hd.find("time", 6) == std::string::npos)
    {
        a = hd;
        b = std::string(initial);
    }
    else
    {
        a = std::string(initial);
        fgets(finish, sizeof(finish), power);
        b = std::string(finish);
    }
    fclose(power);

    // Position of dram_joules
    std::vector<int> start, end;
    // Position of time
    int start_time1, end_time1, start_time2, end_time2;

    // Usage
    double usage = 0;
    double time = 0;

    size_t pos1 = 0;
    size_t pos2 = 0;
    int cnt = 0;

    // Get time
    pos1 = a.find(" ", pos1);
    start_time1 = pos1;
    pos1++;
    pos1 = a.find(" ", pos1);
    end_time1 = pos1;

    pos2 = b.find(" ", pos2);
    start_time2 = pos2;
    pos2++;
    pos2 = b.find(" ", pos2);
    end_time2 = pos2;

    while ((pos2 = b.find(" ", pos2)) != std::string::npos)
    {
        cnt++;

        if (cnt == 4)
        {
            start.push_back(pos2);
        }
        else if (cnt == 5)
        {
            end.push_back(pos2);
            cnt = 0;
        }
        pos2++;
    }

    // Sum usage of all sockets for node level
    if (socket == "-")
    {
        for (size_t i = 0; i < start.size(); i++)
        {
            usage += stod(b.substr(start[i], end[i] - start[i]));
        }
    }
    else
    { // Get individual usage
        int numS = stoi(socket);

        usage = stod(b.substr(start[numS], end[numS] - start[numS]));
    }

    time = stod(b.substr(start_time2, end_time2 - start_time2));
    time = time - stod(a.substr(start_time1, end_time1 - start_time1));
    time = time / 1000; // in second

    std::filesystem::remove(node);

    return std::make_pair(usage, usage / time);
}

std::vector<double> parse_poll_power_idle(const std::string &node)
{

    FILE *power = fopen(node.c_str(), "r");
    fseek(power, 0, SEEK_SET);

    // 25 for powmon and timestamp, 30 for recording per metric
    int buffer_size = 25 + 30 * 5 * 4;
    char header[buffer_size];
    char initial[buffer_size];
    char finish[buffer_size];

    fgets(header, sizeof(header), power);

    fgets(initial, sizeof(initial), power);

    std::string hd = std::string(header);
    std::string a;
    std::string b;

    // Header not printed
    if (hd.find("time", 6) == std::string::npos)
    {
        a = hd;
        b = std::string(initial);
    }
    else
    {
        a = std::string(initial);
        fgets(finish, sizeof(finish), power);
        b = std::string(finish);
    }
    fclose(power);

    // Position of pkg_joules n dram_joules
    std::vector<int> start, end;
    int start_time1, end_time1, start_time2, end_time2;

    size_t pos1 = 0;
    size_t pos2 = 0;
    int cnt = 0;

    // Get time
    pos1 = a.find(" ", pos1);
    start_time1 = pos1;
    pos1++;
    pos1 = a.find(" ", pos1);
    end_time1 = pos1;

    pos2 = b.find(" ", pos2);
    start_time2 = pos2;
    pos2++;
    pos2 = b.find(" ", pos2);
    end_time2 = pos2;

    while ((pos2 = b.find(" ", pos2)) != std::string::npos)
    {
        cnt++;
        if (cnt == 1)
        {
            start.push_back(pos2);
        }
        else if (cnt == 2)
        {
            end.push_back(pos2);
        }
        else if (cnt == 4)
        {
            start.push_back(pos2);
        }
        else if (cnt == 5)
        {
            end.push_back(pos2);
            cnt = 0;
        }
        pos2++;
    }

    // Always in size of 2*#sockets
    std::vector<double> result;
    result.clear();

    double time = stod(b.substr(start_time2, end_time2 - start_time2)) - stod(a.substr(start_time1, end_time1 - start_time1));
    time = time / 1000;

    for (size_t i = 0; i < start.size(); i++)
    {
        double usage = stod(b.substr(start[i], end[i] - start[i]));
        result.push_back(usage / time);
    }

    std::filesystem::remove(node);

    return result;
}

void record_power(const std::string &node, const char *mode)
{
    FILE *in = fopen(node.c_str(), mode);
    fseek(in, 0, SEEK_END);

    if (variorum_poll_power(in))
    {
        std::cerr << "Variorum get power failed!" << std::endl;
        exit(-1);
    }

    fclose(in);
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

    // dram power limit is wrong. so we will cap each socket in the node.
    return socket_power_limit(hostname);
}

double mean(const std::vector<double> &v)
{
    double sum = 0;
    for (double num : v)
    {
        sum += num;
    }
    return sum / v.size();
}

void standard_deviation(const std::vector<double> &data, double &avg, double &sd)
{
    avg = mean(data);

    for (size_t i = 0; i < data.size(); ++i)
        sd += pow(data[i] - avg, 2);

    sd = sqrt(sd / data.size());
}

double covariance(const std::vector<double> &x, const std::vector<double> &y, double mean_x, double mean_y)
{
    double cov = 0;
    for (size_t i = 0; i < x.size(); i++)
    {
        cov += (x[i] - mean_x) * (y[i] - mean_y);
    }
    return cov / x.size();
}

double variance(const std::vector<double> &v, double mean)
{
    double var = 0;
    for (double num : v)
    {
        var += (num - mean) * (num - mean);
    }
    return var / v.size();
}

std::pair<double, double> linear_regression(const std::vector<double> &x, const std::vector<double> &y)
{
    double mean_x = mean(x);
    double mean_y = mean(y);

    double m = covariance(x, y, mean_x, mean_y) / variance(x, mean_x);
    double c = mean_y - m * mean_x;

    return {m, c};
}