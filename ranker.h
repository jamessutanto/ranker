#pragma once

#include <string>
#include <vector>

//int init_topo(const char*);
void rank_dgemm_core(const int, const std::string&);
void rank_dgemm_socket(const int, const std::string&);
void rank_dgemm_node(const int, const std::string&);
void rank_stream_socket(const std::string&);
void rank_stream_node(const std::string&);
void rank_srandom_socket(const std::string&);
void rank_srandom_node(const std::string&);
void rank_core();
void rank_socket();
void rank_node();