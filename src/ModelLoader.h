#include <iostream>
#include <fstream>
#include <map>
#include "nlohmann/json.hpp"



float* loadSchema(const char* path, int* numOfModels, int* modelSize, std::string schemaName);