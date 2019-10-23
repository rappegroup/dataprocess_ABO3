#ifndef interface_h
#define interface_h
#include <iostream>
#include <string>
#include <vector>
std::vector<std::string> split(std::string s,std::string delimiter);
void info(int& cell,std::string& dumpfile,std::string& calistfile,int& velocity_on,int& polarization_on,double& temperature,int& position_variance_on,int& local_die);
#endif
