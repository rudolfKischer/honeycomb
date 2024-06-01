#include "ShaderProgram.h"
#include <glad/glad.h>
#include <fstream>
#include <sstream>
#include <iostream>

const std::string shaderDirPathAbs = "/Users/rudolfkischer/Projects/honeycomb/src/shaders";


ShaderProgram::ShaderProgram(const char* vertexShaderFile, const char* fragmentShaderFile) {
  // build paths
  std::string vertexShaderPath = shaderDirPathAbs + "/" + vertexShaderFile;
  std::string fragmentShaderPath = shaderDirPathAbs + "/" + fragmentShaderFile;
  // build shaders
  unsigned int vertexShader = ShaderProgram::buildShader(vertexShaderPath.c_str(), GL_VERTEX_SHADER);
  unsigned int fragmentShader = ShaderProgram::buildShader(fragmentShaderPath.c_str(), GL_FRAGMENT_SHADER);
  // build program
  ID = ShaderProgram::buildShaderProgram(vertexShader, fragmentShader);

}

ShaderProgram::~ShaderProgram() {
    // glDeleteProgram(ID);
}


unsigned int ShaderProgram::buildShader(const char* shaderPath, GLenum shaderType) {
    std::string shaderSource = ShaderProgram::loadShaderSource(shaderPath);
    const char* shaderSourceCStr = shaderSource.c_str();
    unsigned int shader = glCreateShader(shaderType);
    glShaderSource(shader, 1, &shaderSourceCStr, NULL);
    glCompileShader(shader);
    ShaderProgram::checkShaderCompilation(shader);
    return shader;
}

unsigned int ShaderProgram::buildShaderProgram(unsigned int vertexShader, unsigned int fragmentShader) {
    unsigned int shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);
    ShaderProgram::checkProgramLinking(shaderProgram);
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);
    return shaderProgram;
}

std::string ShaderProgram::loadShaderSource(const char* filePath) {
    std::ifstream shaderFile;
    std::stringstream shaderStream;

    // Ensure ifstream objects can throw exceptions:
    shaderFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
    try {
        // Open file
        shaderFile.open(filePath);
        // Read file's buffer contents into stream
        shaderStream << shaderFile.rdbuf();        
        // Close file handler
        shaderFile.close();
        // Convert stream into string
        return shaderStream.str();
    }
    catch (std::ifstream::failure& e) {
        std::cerr << "ERROR::SHADER::FILE_NOT_SUCCESSFULLY_READ: " << e.what() << '\n';
        return "";
    }
}

void ShaderProgram::checkShaderCompilation(GLuint shader) {
    GLint success;
    GLchar infoLog[512];
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
    if (!success) {
        glGetShaderInfoLog(shader, 512, NULL, infoLog);
        std::cerr << "ERROR::SHADER::COMPILATION_FAILED\n" << infoLog << '\n';
    }
}

void ShaderProgram::checkProgramLinking(GLuint program) {
    GLint success;
    GLchar infoLog[512];
    glGetProgramiv(program, GL_LINK_STATUS, &success);
    if (!success) {
        glGetProgramInfoLog(program, 512, NULL, infoLog);
        std::cerr << "ERROR::SHADER::PROGRAM::LINKING_FAILED\n" << infoLog << '\n';
    }
}

unsigned int makeShaderProgram(const char* vertexShaderFile, const char* fragmentShaderFile) {
    return ShaderProgram(vertexShaderFile, fragmentShaderFile).ID;
}