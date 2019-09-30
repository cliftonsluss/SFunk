#include <fstream>
#include <iostream>
#include <map>

#ifndef READ_CONFIG_H
#define READ_CONFIG_H

class Read_config {
  public:
    Read_config(std::string &config_file,
                std::map<std::string, std::string> &config_map);
      void get_config();

  private:
    std::ifstream config_file;
    std::map<std::string, std::string> config_map;



};


#endif
