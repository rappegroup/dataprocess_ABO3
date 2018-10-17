#ifndef interface_h
#define interface_h
#include <iostream>
#include <string>
#include <vector>
std::vector<std::string> split(std::string s,std::string delimiter);
void info(int& cell,std::string& dumpfile,std::string& calistfile,bool& velocity_on,bool& polarization_on);
#endif
