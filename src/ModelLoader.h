




float* loadModel(const char* path, int* numOfModels);

struct Sphere {
    float position[3];
    float radius;
    float color[3];
    float emission[3];
    float orm[3];
    // o = opacity
    // r = roughness
    // m = metallic
};

//     // Pack spheres into a float array
void packSpheres(Sphere* spheres, int numSpheres, float* data);