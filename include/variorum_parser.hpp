#pragma once

#include <iostream>
#include <vector>

double parse_poll_power_pkg(const std::string &node, const std::string &socket);
double parse_poll_power_dram(const std::string &node, const std::string &socket);
std::vector<double> parse_poll_power_idle(const std::string &node);
void record_power(const std::string &node, const char *mode);
int socket_power_limit(const char *node);
int node_power_limit(const char *node);