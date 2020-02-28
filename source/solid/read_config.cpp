#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>
#include "read_config.h"

Read_config::Read_config(std::string &config_file){
       // std::map<std::string,std::string> &config_map) {
  // Read_config::config_file.open(config_file, &std::ifstream::in);
  Read_config::config_file = config_file;
  // Read_config::config_map = config_map;
}


void Read_config::get_config(){
  // if(config_file.is_open())
  // {
    std::ifstream i(Read_config::config_file);
    i >> Read_config::j;
  // }


  // std::ifstream is RAII, i.e. no need to call close
  //   if(config_file.is_open())
  //   {
  //     std::string line;
  //     while(getline(config_file, line)){
  //       line.erase(std::remove_if(line.begin(), line.end(), isspace),
  //                             line.end());
  //       if(line[0] == '#' || line.empty())
  //           continue;
  //       auto delimiterPos = line.find("=");
  //       config_map[line.substr(0, delimiterPos)] = line.substr(delimiterPos + 1);
  //     }
  //   }
  //   else {
  //     std::cerr << "Couldn't open config file for reading.\n";
  //   }
}
