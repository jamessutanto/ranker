#pragma once

#include <iostream>
#include <vector>

struct Row
{
    std::vector<std::string> columns;
};

std::istream &operator>>(std::istream &str, Row &data);
bool compare_last_reversed(const Row &r1, const Row &r2);
bool compare_secondlast(const Row &r1, const Row &r2);
void to_vector_header(const std::string &filename, std::vector<Row> &local);
void to_vector(const std::string &filename, std::vector<Row> &local);
void write_stream_result(const std::string &path_to, const std::string &path_from, const std::string &node, const std::string &socket, const std::string &threads, const int &cache_line_size);
void write_nas_result(const std::string &path_to, const std::string &path_from, const std::string &node, const std::string &socket, const std::string &threads, const int &cache_line_size);
void write_nas_core_result(const std::string &path_to, const std::string &path_from, const std::string &node, const std::string &socket, const std::string &threads);
void calculate_score(std::vector<Row> &result_e, std::vector<std::vector<Row>> &result_t, const std::vector<Row> &data, const std::vector<Row> &data_nolimit, const std::vector<Row> &nas_data, const std::string &benchmark);
void calculate_efficiency(std::vector<Row> &result, const std::vector<Row> &data, const std::vector<Row> &data_nolimit, const std::vector<Row> &ep_data, const std::string &benchmark);
void calculate_dram_dynamic(const std::string &path, const std::string &final);

void change_affinity(const std::vector<int> &cpu);
void set_omp(const std::vector<int> &cpu);

std::pair<double, double> parse_poll_power_pkg(const std::string &node, const std::string &socket);
std::pair<double, double> parse_poll_power_dram(const std::string &node, const std::string &socket);
std::vector<double> parse_poll_power_idle(const std::string &node);
void record_power(const std::string &node, const char *mode);
int get_socket_power_limit();
int get_node_power_limit();
int socket_power_limit(const char *node);
int node_power_limit(const char *node);
int get_n_socket();

void standard_deviation(const std::vector<double> &data, double &avg, double &sd);
std::pair<double, double> linear_regression(const std::vector<double>& x, const std::vector<double>& y);