#include "ModelLoader.h"
#include "OBJ_Loader.h"
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

std::vector<Entry> MeshEntries = {
    {"position", 3},
    {"material", 1},
    {"obj_path_str", 0}
};


Schema sphereSchema = initSchema(sphereEntries, "spheres");
Schema materialSchema = initSchema(materialEntries, "materials");
Schema aabbSchema = initSchema(aabbEntries, "aabb");
Schema meshSchema = initSchema(MeshEntries, "meshes");


// map from string to schema
std::map<std::string, Schema> schemas = {
    {"spheres", sphereSchema},
    {"materials", materialSchema},
    {"aabb", aabbSchema},
    {"meshes", meshSchema}
};

// we want to be able to load mesh objects
// mesh objects will contain the normal defiintions relative to the scene, like the position of the object, the material
// and maybe rotation and the like
// they will also have the path to the obj file
// we need to be able to load the the obj file into a vbo texture
// we then need to return a list of mesh objects
// the mesh objects do not contains the verts, or idices, they just contain where they start and how many there are
float* loadMeshes(const char* path, int* num_meshes, int* mesh_size, int* num_verts, int* num_faces) {
    Schema schema = schemas["meshes"];
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

    // the way this will work is that the data return back will be as such
    // | mesh1 | mesh2 ... | VERTICES | FACE INDICES |
    // we will first put vertices and indices into a seperate array
    // and then append them to the end of the mesh array

    objl::Loader obj_loader;
    std::vector<float> vertices;
    std::vector<int> indices;

    int vert_offset = 0;
    int index_offset = 0;

    int cur_num_meshes = 0;


    std::vector<float> data;
    for (const auto& item : j[schema.name]) {
        // we first need to get the path to the obj file
        std::string obj_path = item["obj_path_str"];
        // we then need to load the obj file into a vbo
        bool loadout = obj_loader.LoadFile(obj_path);
        if (!loadout) {
            std::cerr << "Unable to load obj file: " << obj_path << std::endl;
            return nullptr;
        }

        for (int i =0; i < obj_loader.LoadedMeshes.size(); i++) {
            objl::Mesh cur_mesh = obj_loader.LoadedMeshes[i];

            // std::cout << "Mesh: " << i << " " << cur_mesh.MeshName << std::endl;

            // debug current mesh
            // for (int j = 0; j < cur_mesh.Vertices.size(); j++)
            // {
            //   std::cout << "V" << j << ": " <<
            //     "P(" << cur_mesh.Vertices[j].Position.X << ", " << cur_mesh.Vertices[j].Position.Y << ", " << cur_mesh.Vertices[j].Position.Z << ") " <<
            //     "N(" << cur_mesh.Vertices[j].Normal.X << ", " << cur_mesh.Vertices[j].Normal.Y << ", " << cur_mesh.Vertices[j].Normal.Z << ") " <<
            //     "TC(" << cur_mesh.Vertices[j].TextureCoordinate.X << ", " << cur_mesh.Vertices[j].TextureCoordinate.Y << ")\n";
            // }



            // we will now load the mesh data
            // we want to get the total number of vertices and faces
            int num_verts = cur_mesh.Vertices.size();
            int num_faces = cur_mesh.Indices.size();

            // MESH STRUCTURE
            // | position | material | vert_offset | index_offset | num_verts | num_faces | minp | maxp 
            // we will now load the mesh data
            // MESH SIZE = 3 + 1 + 1 + 1 + 1 + 1 + 3 + 3 = 14

            // we first load all the schema entries from the json, except for the obj_path_str
            for (int i = 0; i < schema.fields.size(); i++) {
                if (schema.fields[i].name == "obj_path_str") {
                    continue;
                }
                for (int k = 0; k < schema.fields[i].size; k++) {
                    if (schema.fields[i].size == 1) {
                        data.push_back(item[schema.fields[i].name]);
                    } else {
                        data.push_back(item[schema.fields[i].name][k]);
                    }
                }
            }

            // then we add the vert_offset, index_offset, num_verts, and num_faces
            data.push_back(vert_offset);
            data.push_back(index_offset);
            data.push_back(num_verts);
            data.push_back(num_faces);

            // print the number of vertices
            std::cout << "Num Verts: " << num_verts << std::endl;
            std::cout << "Num Faces: " << num_faces << std::endl;






            // then we load the vertices and indices into the arrays
            for (int i = 0; i < num_verts; i++) {
                vertices.push_back(cur_mesh.Vertices[i].Position.X);
                vertices.push_back(cur_mesh.Vertices[i].Position.Y);
                vertices.push_back(cur_mesh.Vertices[i].Position.Z);
            }

            // go through the vertices
            float minp[3] = {9999999, 9999999, 9999999};
            float maxp[3] = {-9999999, -9999999, -9999999};
            for (int i = 0; i < num_verts; i++) {
              for (int j = 0; j < 3; j++) {
                if (vertices[i * 3 + j] < minp[j]) {
                  minp[j] = vertices[i * 3 + j];
                }
                if (vertices[i * 3 + j] > maxp[j]) {
                  maxp[j] = vertices[i * 3 + j];
                }
              }
            }

            for (int i = 0; i < 3; i++) {
              data.push_back(minp[i]);
            }
            for (int i = 0; i < 3; i++) {
              data.push_back(maxp[i]);
            }

            // print min and max p
            std::cout << "Minp: " << minp[0] << " " << minp[1] << " " << minp[2] << std::endl;
            std::cout << "Maxp: " << maxp[0] << " " << maxp[1] << " " << maxp[2] << std::endl;


            for (int i = 0; i < num_faces; i++) {
                indices.push_back(cur_mesh.Indices[i]);
            }

            // we then update the vert_offset and index_offset
            vert_offset += num_verts;
            index_offset += num_faces;

            // we now update the number of meshes loaded
            cur_num_meshes++;


        }


    }

    // we now need to append the vertices and indices to the end of the data array
    for (int i = 0; i < vertices.size(); i++) {
        data.push_back(vertices[i]);
    }

    // we now need to append the indices to the end of the data array
    for (int i = 0; i < indices.size(); i++) {
        data.push_back(indices[i]);
    }

    *num_verts = vertices.size() / 3;
    *num_faces = indices.size() / 3;
    *num_meshes = cur_num_meshes;
    *mesh_size = schema.size + 4 + 6;
    float* dataArr = new float[data.size()];
    for (int i = 0; i < data.size(); i++) {
        dataArr[i] = data[i];
    }
    return dataArr;
}








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