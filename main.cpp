
#include "ranker.h"

#include <vector>
#include <string>

int main()
{
    std::vector<std::string> path = {"/u/home/sutantoj/ba/ranker/topology.xml"};
    std::string result = "/u/home/sutantoj/ba/ranker/result/res.csv";
    int powerLimit = 100;
    rank_dgemm_socket(powerLimit, result);

    return 0;
}