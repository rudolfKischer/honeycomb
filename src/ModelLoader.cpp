#include "ModelLoader.h"
#include <iostream>
#include <fstream>
#include "nlohmann/json.hpp" // nlohmann

using json = nlohmann::json;

const int datum = 13;

void packSpheres(Sphere* spheres, int numSpheres, float* data) {
    
    for (int i = 0; i < numSpheres; i++) {
        data[i * datum] = spheres[i].position[0];
        data[i * datum + 1] = spheres[i].position[1];
        data[i * datum + 2] = spheres[i].position[2];
        data[i * datum + 3] = spheres[i].radius;
        data[i * datum + 4] = spheres[i].color[0];
        data[i * datum + 5] = spheres[i].color[1];
        data[i * datum + 6] = spheres[i].color[2];
        data[i * datum + 7] = spheres[i].emission[0];
        data[i * datum + 8] = spheres[i].emission[1];
        data[i * datum + 9] = spheres[i].emission[2];// orm
        data[i * datum + 10] = spheres[i].orm[0];
        data[i * datum + 11] = spheres[i].orm[1];
        data[i * datum + 12] = spheres[i].orm[2];
    }
}

float* loadModel(const char* path, int* numOfModels) {

    // load the model
        // Read the JSON file
    std::ifstream file(path);
    if (!file) {
        std::cerr << "Unable to open file: " << path << std::endl;
        return nullptr;
    }
    
    json j;
    file >> j;
    file.close();
    
    // Check if the JSON contains "spheres"
    if (!j.contains("spheres") || !j["spheres"].is_array()) {
        std::cerr << "JSON does not contain a valid 'spheres' array." << std::endl;
        return nullptr;
    }
    
    // Extract spheres
    std::vector<Sphere> spheres;
    for (const auto& item : j["spheres"]) {
        Sphere sphere;
        item.at("position").get_to(sphere.position);
        item.at("radius").get_to(sphere.radius);
        item.at("color").get_to(sphere.color);
        item.at("emission").get_to(sphere.emission);
        item.at("orm").get_to(sphere.orm);
        spheres.push_back(sphere);
    }

    *numOfModels = spheres.size();
    
    // Allocate memory for the packed data
    float* data = new float[spheres.size() * datum];
    packSpheres(spheres.data(), spheres.size(), data);
    
    return data;

}