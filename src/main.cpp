
#include "ShaderProgram.h"
#include "WindowManager.h"
#include "RenderTexture.h"
#include "Quad.h"
#include "ModelLoader.h"

#include <iostream>
#include <random>
#include <string>
#include <chrono>

std::random_device rd;  // obtain a random number from hardware
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis(0.0, 1.0);

const unsigned int TEX_WIDTH = 1000;
const unsigned int TEX_HEIGHT = 1000;
int SCR_WIDTH = 1000;
int SCR_HEIGHT = 1000;

const int datum = 13;
float camPos[3] = {0.0, 1.5, 3.0};
float speed = 0.01;

void setRenderTarget(unsigned int shaderProgram, unsigned int framebuffer, int w, int h);
void updateCameraPosition(GLFWwindow* window);
void fpsCounter(double* lastTime, int* nbFrames);
bool getCamMoved(float* camPos, float* prevCamPos);
void setShaderProgramUniforms(unsigned int shaderProgram, int frameCount, int TEX_WIDTH, int TEX_HEIGHT, int numOfModels, float* sphereData, float* camPos, bool camMoved);



int main()
{

    //load model test
    std::string assetsFolder = "/Users/rudolfkischer/Projects/honeycomb/assets/";
    std::string testModel = "test.json";
    std::string cornellBox = "cornell.json";
    std::string cornellBoxTwoSpheres = "cornellTwoSpheres.json";
    std::string grid = "grid.json";


    // std::string modelPathStr = assetsFolder + testModel;
    // std::string modelPathStr = assetsFolder + cornellBox;
    // std::string modelPathStr = assetsFolder + grid;
    std::string modelPathStr = assetsFolder + cornellBoxTwoSpheres;


    const char* modelPath = (const char*)modelPathStr.c_str();
    int numOfModels;
    float* sphereData = loadModel(modelPath, &numOfModels);
    if (sphereData == nullptr) {
        std::cout << "Failed to load model" << std::endl;
        return -1;
    }

    // Initialize window
    GLFWwindow* window = make_window(SCR_WIDTH, SCR_HEIGHT, "LearnOpenGL");
    
    // Initialize GLAD
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    // Load shaders
    unsigned int shaderProgram = makeShaderProgram("vertex.glsl", "fragment.glsl");
    unsigned int screenShaderProgram = makeShaderProgram("screenVertex.glsl", "screenFragment.glsl");
    std::cout << "OpenGL Version: " << glGetString(GL_VERSION) << std::endl;

    
    Quad* quad = new Quad(); // quad is used to put texture on
    
    RenderTexture* renderTexture = new RenderTexture(TEX_WIDTH, TEX_HEIGHT); //this holds both a frame buffer, and a texture the frame buffer is what is rendered to, the texture object is what stores the render

    // Time control for pixel lighting
    double lastTime = glfwGetTime();
    int nbFrames = 0;
    int frameCount = 0;

    float prevCamPos[3];
    for (int i = 0; i < 3; i++) {
        prevCamPos[i] = camPos[i];
    }
    

    // Render loop
    while (!glfwWindowShouldClose(window))
    {    

        updateCameraPosition(window);
        glfwGetFramebufferSize(window, &SCR_WIDTH, &SCR_HEIGHT); // update screen size

        bool camMoved = getCamMoved(camPos, prevCamPos);
        if (camMoved) { frameCount = 0; } else { frameCount++; }
        fpsCounter(&lastTime, &nbFrames);

        setRenderTarget(shaderProgram, renderTexture->FBO, TEX_WIDTH, TEX_HEIGHT); // render to canvas
        setShaderProgramUniforms(shaderProgram, frameCount, TEX_WIDTH, TEX_HEIGHT, numOfModels, sphereData, camPos, camMoved); // set shader uniforms

        glActiveTexture(GL_TEXTURE1); // bind previous frame texture
        glBindTexture(GL_TEXTURE_2D, renderTexture->prevTexture);
        glUniform1i(glGetUniformLocation(shaderProgram, "prevFrame"), 1);
        quad->draw();
        renderTexture->copy(); // copy the rendered texture to the previous frame texture

        // Render to screen
        quad->texture = renderTexture->texture;
        setRenderTarget(screenShaderProgram, 0, SCR_WIDTH, SCR_HEIGHT);
        quad->draw();

        // apply changes to window and get input events
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // cleanup shaders
    glDeleteProgram(shaderProgram);
    glDeleteProgram(screenShaderProgram);
    glfwTerminate();
    return 0;
}

void setRenderTarget(unsigned int shaderProgram, unsigned int framebuffer, int w, int h) {
    // Bind to framebuffer and draw scene to texture
    glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);
    glViewport(0, 0, w, h);
    glClear(GL_COLOR_BUFFER_BIT);
    glUseProgram(shaderProgram);
}

//keyboard callback
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{ 
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
}

void updateCameraPosition(GLFWwindow* window) {

    // if control is pressed, up the
    if (glfwGetKey(window, GLFW_KEY_LEFT_CONTROL) == GLFW_PRESS)
        speed = 0.1;
    else
        speed = 0.01;

    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        camPos[2] -= speed;
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        camPos[2] += speed;
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        camPos[0] -= speed;
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        camPos[0] += speed;
    if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS)
        camPos[1] += speed;
    if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS)
        camPos[1] -= speed;
}

void fpsCounter(double* lastTime, int* nbFrames) {
    double currentTime = glfwGetTime();
    (*nbFrames)++;
    if (currentTime - *lastTime >= 1.0) { // If last print was more than 1 sec ago
        std::cout << *nbFrames << " FPS" << std::endl;
        *nbFrames = 0;
        *lastTime += 1.0;
    }
}

bool getCamMoved(float* camPos, float* prevCamPos) {
    bool camMoved = false;
    for (int i = 0; i < 3; i++) {
        if (abs(camPos[i] - prevCamPos[i]) > 0.0001) 
            camMoved = true;
        prevCamPos[i] = camPos[i];
    }
    return camMoved;
}

void setShaderProgramUniforms(unsigned int shaderProgram, int frameCount, int TEX_WIDTH, int TEX_HEIGHT, int numOfModels, float* sphereData, float* camPos, bool camMoved) {
    glUniform1i(glGetUniformLocation(shaderProgram, "timePassed"), frameCount);
    glUniform1i(glGetUniformLocation(shaderProgram, "texWidth"), TEX_WIDTH);
    glUniform1i(glGetUniformLocation(shaderProgram, "texHeight"), TEX_HEIGHT);
    glUniform1i(glGetUniformLocation(shaderProgram, "numSpheres"), numOfModels);
    glUniform1fv(glGetUniformLocation(shaderProgram, "spheredata"), numOfModels * datum, sphereData);
    glUniform3fv(glGetUniformLocation(shaderProgram, "camPos"), 1, camPos);
    glUniform1i(glGetUniformLocation(shaderProgram, "camMoved"), camMoved);
    glUniform1f(glGetUniformLocation(shaderProgram, "seed1"), dis(gen));
    glUniform1f(glGetUniformLocation(shaderProgram, "uTime"), glfwGetTime());
}



