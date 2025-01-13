#pragma once

#include <iostream>

struct AppendEndL {
  ~AppendEndL() { 
    std::cerr << std::endl;  
    if ( exit_code ) {
      exit(exit_code);
    }
  }
  AppendEndL(int c) {
    exit_code = c;
  }
  int exit_code = 0;
};

struct AppendEndLAndExit {
  ~AppendEndLAndExit() { 
    std::cerr << std::endl; 
    exit(-1); 
  }
};

/*
#define INFO "I"
#define WARNING "W"
#define ERROR "E"
*/

#define LOG(s) (AppendEndL(std::string(#s) != "ERROR" ? 0 : -1), std::cerr << __FILE__ << ":" << __LINE__ << " " << #s << " ")

#define CHECK(c) if ( !(c) ) (AppendEndLAndExit(), std::cerr << __FILE__ << ":" << __LINE__ << " " << #c << " failed ")

#define TF_CHECK_OK(expr) CHECK(expr.ok())

#define CHECK_GE(a, b) CHECK((a) >= (b))
#define CHECK_GT(a, b) CHECK((a) > (b))
#define CHECK_LE(a, b) CHECK((a) <= (b))
