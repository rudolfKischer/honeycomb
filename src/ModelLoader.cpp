#include "ModelLoader.h"
 // nlohmann

using json = nlohmann::json;

struct Entry {
    std::string name;
    int size;
};

// schemas will be mappings from the struct to a flattened array
struct Field {
    std::string name;
    int size;
    int offset;
};
// we will use the schema to flatten the struct into a float array

struct Schema {
    std::string name;
    std::vector<Field> fields;
    int size;
};


// offset, and the total size of the schema should be calculated
Schema initSchema(std::vector<Entry> entries, std::string name) {
    Schema schema;
    schema.name = name;
    schema.size = 0;
    for (int i = 0; i < entries.size(); i++) {
        Field field;
        field.name = entries[i].name;
        field.size = entries[i].size;
        field.offset = schema.size;
        schema.size += entries[i].size;
        schema.fields.push_back(field);
    }
    return schema;
}

std::vector<Entry> sphereEntries = {
    {"position", 3},
    {"radius", 1},
    // {"color", 3},
    // {"emission", 3},
    // {"orm", 3},
    {"material", 1}
};

std::vector<Entry> materialEntries = {
    {"color", 3},
    {"emission", 3},
    {"orm", 3}
};


std::vector<Entry> aabbEntries = {
    {"minp", 3},
    {"maxp", 3},
    {"material", 1}
};


Schema sphereSchema = initSchema(sphereEntries, "spheres");
Schema materialSchema = initSchema(materialEntries, "materials");
Schema aabbSchema = initSchema(aabbEntries, "aabb");


// map from string to schema
std::map<std::string, Schema> schemas = {
    {"spheres", sphereSchema},
    {"materials", materialSchema},
    {"aabb", aabbSchema}
};


float* loadSchema(const char* path, int* num, int* size, std::string schemaName) {
    Schema schema = schemas[schemaName];
    // Read the JSON file
  std::ifstream file(path);
    if (!file) {
        std::cerr << "Unable to open file: " << path << std::endl;
        return nullptr;
    }
    
    json j;
    file >> j;
    file.close();
    
    // Check if the JSON contains Schema.name
    if (!j.contains(schema.name) || !j[schema.name].is_array()) {
        std::cerr << "JSON does not contain a valid '" << schema.name << "' array." << std::endl;
        return nullptr;
    }

    std::vector<float> data;
    for (const auto& item : j[schema.name]) {
        for (int i = 0; i < schema.fields.size(); i++) {
            for (int k = 0; k < schema.fields[i].size; k++) {
                if (schema.fields[i].size == 1) {
                    // print out the name of the current field and the value
                    // std::cout << schema.fields[i].name << ": " << item[schema.fields[i].name] << std::endl;
                    data.push_back(item[schema.fields[i].name]);
                } else {
                    data.push_back(item[schema.fields[i].name][k]);
                }
            }
        }
    }

    *num = j[schema.name].size();
    float* dataArr = new float[data.size()];
    for (int i = 0; i < data.size(); i++) {
        dataArr[i] = data[i];
    }
    *size = schema.size;
    return dataArr;
}