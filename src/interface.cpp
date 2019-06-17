#include <string>
#include <iostream>
#include <vector>
#include <sstream>
std::vector<std::string> split(std::string s, std::string delimiter) {
	    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
			std::string token;
			std::vector<std::string> res;
			while ((pos_end = s.find (delimiter, pos_start)) != std::string::npos) {
			 token = s.substr(pos_start, pos_end - pos_start);
			 pos_start = pos_end + delim_len;
			 res.push_back(token);
			}
			res.push_back(s.substr(pos_start));
		  return res;
}
void info(int& cell,std::string& dumpfile,std::string& calistfile,bool& velocity_on,bool& polarization_on,double& temperature){
	std::string temp;
	getline(std::cin,temp);//rm the first line, this is the tile and useless.
	std::vector<std::string> input;
	std::string flag;
	std::istringstream stre;
	getline(std::cin,temp);
	while(temp.find("/")==std::string::npos){
		input=split(temp,"=");
		flag=input[1];
		if(input[0].find("filename")!=std::string::npos){
			flag=flag.substr(0,flag.find(","));
			stre.str(flag);
		  stre>>dumpfile;
			stre.clear();
		}
		else if(input[0].find("cell")!=std::string::npos){
			flag=flag.substr(0,flag.find(","));
			stre.str(flag);
			stre>>cell;
			stre.clear();
		}
		else if(input[0].find("solution_list")!=std::string::npos){
			flag=flag.substr(0,flag.find(","));
			stre.str(flag);
			stre>>calistfile;
			stre.clear();
		}
		else if(input[0].find("auto_velocity")!=std::string::npos){
			flag=flag.substr(0,flag.find(","));
			stre.str(flag);
			stre>>velocity_on;
			stre.clear();
		}
		else if(input[0].find("polarization")!=std::string::npos){
			flag=flag.substr(0,flag.find(","));
			stre.str(flag);
			stre>>polarization_on;
			stre.clear();
		}
		else if(input[0].find("temp")!=std::string::npos){
			flag=flag.substr(0,flag.find(","));
			stre.str(flag);
			stre>>temperature;
			stre.clear();
		}
	getline(std::cin,temp);
	}
}
