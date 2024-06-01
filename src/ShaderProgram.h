#pragma once
#include <glad/glad.h>
#include <string>

class ShaderProgram {
public:
    unsigned int ID;
    ShaderProgram(const char* vertexShaderFile, const char* fragementShaderFile);
    ~ShaderProgram();

private:
    static std::string loadShaderSource(const char* filePath);
    static void checkShaderCompilation(GLuint shader);
    static void checkProgramLinking(GLuint program);
    static unsigned int buildShader(const char* shaderPath, GLenum shaderType);
    static unsigned int buildShaderProgram(unsigned int vertexShader, unsigned int fragmentShader);
};

unsigned int makeShaderProgram(const char* vertexShaderFile, const char* fragmentShaderFile);

