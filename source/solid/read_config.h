#include <fstream>
#include <iostream>
#include <map>
#include "json.h"

#ifndef READ_CONFIG_H
#define READ_CONFIG_H

using json = nlohmann::json;

class Read_config {
  public:
    Read_config(std::string &config_file);
                // std::map<std::string, std::string> &config_map);
    void get_config();
    json j;

  private:
    // std::ifstream config_file;
    // std::map<std::string, std::string> config_map;
    std::string config_file;

  // public:
  //   json j;

};




#endif
