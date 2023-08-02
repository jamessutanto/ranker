#include "variorum_parser.hpp"

#include <vector>
#include <unistd.h>
#include <filesystem>
#include <mpi.h>

extern "C"
{
#include "variorum.h"
}

int socket_power_limit(const char *node)
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

    int start = 0;
    int end = -1;

    size_t pos = 0;
    int cnt = 0;

    while ((pos = a.find(" ", pos)) != std::string::npos)
    {
        cnt++;
        if ((cnt - 3) % 5 == 0)
        {
            start = pos;
        }
        else if ((cnt - 4) % 5 == 0)
        {
            end = pos;
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

    int cnt = 0;
    size_t pos = 0;

    while ((pos = a.find(" ", pos)) != std::string::npos)
    {
        cnt++;
        if ((cnt - 3) % 5 == 0)
        {
            start.push_back(pos);
        }
        else if ((cnt - 4) % 5 == 0)
        {
            end.push_back(pos);
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

double parse_poll_power_pkg(const std::string &node, const std::string &socket)
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
    std::vector<int> start1, end1, start2, end2;

    // Usage
    double before = 0;
    double after = 0;

    size_t pos1 = 0;
    size_t pos2 = 0;
    int cnt = 0;

    // Removing _POWMON and time
    pos1 = a.find(" ", pos1);
    pos1++;
    pos1 = a.find(" ", pos1);

    pos2 = b.find(" ", pos2);
    pos2++;
    pos2 = b.find(" ", pos2);

    while ((pos1 = a.find(" ", pos1)) != std::string::npos)
    {
        pos2 = b.find(" ", pos2);
        cnt++;
        if (cnt == 1)
        {
            start1.push_back(pos1);
            start2.push_back(pos2);
        }
        else if (cnt == 2)
        {
            end1.push_back(pos1);
            end2.push_back(pos2);
        }
        else if (cnt == 5)
        {
            cnt = 0;
        }
        pos1++;
        pos2++;
    }

    std::cout << a << std::endl;
    std::cout << b << std::endl;

    // Sum usage of all sockets for node level
    if (socket == "-")
    {
        for (size_t i = 0; i < start1.size(); i++)
        {
            before += stod(a.substr(start1[i], end1[i] - start1[i]));
            after += stod(b.substr(start2[i], end2[i] - start2[i]));
        }
    }
    else
    { // Get individual usage
        int numS = stoi(socket);

        before = stod(a.substr(start1[numS], end1[numS] - start1[numS]));
        after = stod(b.substr(start2[numS], end2[numS] - start2[numS]));
    }

    std::filesystem::remove(node);

    return after - before;
}

double parse_poll_power_dram(const std::string &node, const std::string &socket)
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
    std::vector<int> start1, end1, start2, end2;
    int start_time1, end_time1, start_time2, end_time2;

    // Usage
    double before = 0;
    double after = 0;

    size_t pos1 = 0;
    size_t pos2 = 0;
    int cnt = 0;

    // Get time
    pos1 = a.find(" ", pos1);
    start_time1 = pos1;
    pos1++;
    pos1 = a.find(" ", pos1);
    end_time1 = pos1;
    pos1++;

    pos2 = b.find(" ", pos2);
    start_time2 = pos2;
    pos2++;
    pos2 = b.find(" ", pos2);
    end_time2 = pos2;
    pos2++;

    while ((pos1 = a.find(" ", pos1)) != std::string::npos)
    {
        pos2 = b.find(" ", pos2);
        cnt++;

        if (cnt == 4)
        {
            start1.push_back(pos1);
            start2.push_back(pos2);
        }
        else if (cnt == 5)
        {
            end1.push_back(pos1);
            end2.push_back(pos2);
            cnt = 0;
        }
        pos1++;
        pos2++;
    }

    std::cout << a << std::endl;
    std::cout << b << std::endl;
    
    double time = stod(a.substr(start_time1, end_time1 - start_time1)) - stod(b.substr(start_time2, end_time2 - start_time2));
    time = time / 1000;

    // Sum usage of all sockets for node level
    if (socket == "-")
    {
        for (size_t i = 0; i < start1.size(); i++)
        {
            before += stod(a.substr(start1[i], end1[i] - start1[i]));
            after += stod(b.substr(start2[i], end2[i] - start2[i]));
        }
    }
    else
    { // Get individual usage
        int numS = stoi(socket);

        before = stod(a.substr(start1[numS], end1[numS] - start1[numS]));
        after = stod(b.substr(start2[numS], end2[numS] - start2[numS]));
    }

    double joule = after - before;

    std::filesystem::remove(node);

    return joule / time;
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
    std::vector<int> start1, end1, start2, end2;
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

    while ((pos1 = a.find(" ", pos1)) != std::string::npos)
    {
        pos2 = b.find(" ", pos2);
        cnt++;
        if (cnt == 1)
        {
            start1.push_back(pos1);
            start2.push_back(pos2);
        }
        else if (cnt == 2)
        {
            end1.push_back(pos1);
            end2.push_back(pos2);
        }
        else if (cnt == 4)
        {
            start1.push_back(pos1);
            start2.push_back(pos2);
        }
        else if (cnt == 5)
        {
            end1.push_back(pos1);
            end2.push_back(pos2);
            cnt = 0;
        }
        pos1++;
        pos2++;
    }

    // Always in size of 2*#sockets
    std::vector<double> result;
    result.clear();

    double time = stod(b.substr(start_time2, end_time2 - start_time2)) - stod(a.substr(start_time1, end_time1 - start_time1));
    time = time / 1000;

    std::cout << a << std::endl;
    std::cout << b << std::endl;

    for (size_t i = 0; i < start1.size(); i++)
    {
        double before = stod(a.substr(start1[i], end1[i] - start1[i]));
        double after = stod(b.substr(start2[i], end2[i] - start2[i]));
        result.push_back((after - before) / time);
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