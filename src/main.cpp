
#include "ShaderProgram.h"
#include "WindowManager.h"
#include "RenderTexture.h"
#include "Quad.h"
#include "ModelLoader.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

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

float camPos[3] = {0.0, 2.5, 6.0};
float camFront[3] = {0.0, 0.0, -1.0};
float speed = 0.01;

void setRenderTarget(unsigned int shaderProgram, unsigned int framebuffer, int w, int h);
void updateCameraPosition(GLFWwindow* window, double deltaTime);
void fpsCounter(double* lastTime, int* nbFrames);
bool getCamMoved(float* camPos, float* prevCamPos, float* camFront, float* prevCamFront);
void setShaderProgramUniforms(unsigned int shaderProgram, int frameCount, int TEX_WIDTH, int TEX_HEIGHT, int numOfModels, float* sphereData, float* camPos, bool camMoved, int modelSize);


GLuint loadTexture(const char* path, int* width, int* height, int* nrChannels) {
    // Check file extension to determine the format
    std::string filePath(path);
    bool isHDR = filePath.substr(filePath.find_last_of(".") + 1) == "hdr";

    unsigned char* data = nullptr;
    float* hdrData = nullptr;

    if (isHDR) {
        hdrData = stbi_loadf(path, width, height, nrChannels, 0);
        if (hdrData) {
            std::cout << "Loaded HDR texture: " << path << std::endl;
        } else {
            std::cout << "Failed to load HDR texture" << std::endl;
            return 0;
        }
    } else {
        data = stbi_load(path, width, height, nrChannels, 0);
        if (data) {
            std::cout << "Loaded texture: " << path << std::endl;
        } else {
            std::cout << "Failed to load texture" << std::endl;
            return 0;
        }
    }

    // TEXTURE LOAD
    GLuint bgTexture;
    glGenTextures(1, &bgTexture);
    glBindTexture(GL_TEXTURE_2D, bgTexture);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    // Set texture filtering parameters
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    if (isHDR) {
        if (*nrChannels == 3)
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB16F, *width, *height, 0, GL_RGB, GL_FLOAT, hdrData);
        else if (*nrChannels == 4)
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, *width, *height, 0, GL_RGBA, GL_FLOAT, hdrData);
        else {
            std::cout << "Failed to load HDR texture" << std::endl;
            stbi_image_free(hdrData);
            return 0;
        }
        stbi_image_free(hdrData);
    } else {
        if (*nrChannels == 3)
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, *width, *height, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
        else if (*nrChannels == 4)
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, *width, *height, 0, GL_RGBA, GL_UNSIGNED_BYTE, data);
        else {
            std::cout << "Failed to load texture" << std::endl;
            stbi_image_free(data);
            return 0;
        }
        stbi_image_free(data);
    }

    glGenerateMipmap(GL_TEXTURE_2D);

    return bgTexture;
}


int main()
{

    //load model test
    std::string assetsFolder = "/Users/rudolfkischer/Projects/honeycomb/assets/models/";
    std::string testModel = "test.json";
    std::string cornellBox = "cornell.json";
    std::string cornellBoxTwoSpheres = "cornellTwoSpheres.json";
    std::string grid = "grid.json";
    std::string simple = "simple.json";
    std::string cornellSunlight = "cornellSunlight.json";


    // std::string modelPathStr = assetsFolder + testModel;
    // std::string modelPathStr = assetsFolder + cornellBox;
    // std::string modelPathStr = assetsFolder + grid;
    std::string modelPathStr = assetsFolder + cornellBoxTwoSpheres;
    // std::string modelPathStr = assetsFolder + cornellSunlight;
    // std::string modelPathStr = assetsFolder + simple;


    const char* modelPath = (const char*)modelPathStr.c_str();
    int numOfModels, modelSize;
    float* sphereData = loadSchema(modelPath, &numOfModels, &modelSize, "spheres");
    if (sphereData == nullptr) {
        std::cout << "Failed to load model" << std::endl;
        return -1;
    }
    int numOfMaterials, materialSize;
    float* materialData = loadSchema(modelPath, &numOfMaterials, &materialSize, "materials");
    if (materialData == nullptr) {
        std::cout << "Failed to load model" << std::endl;
        return -1;
    }
    int numOfAABBs, aabbSize;
    float* aabbData = loadSchema(modelPath, &numOfAABBs, &aabbSize, "aabb");
    if (aabbData == nullptr) {
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
    float prevCamFront[3];
    for (int i = 0; i < 3; i++) {
        prevCamPos[i] = camPos[i];
        prevCamFront[i] = camFront[i];
    }


    int bgWidth, bgHeight, bgChannels;
    // std::string bgPath = "/Users/rudolfkischer/Projects/honeycomb/assets/textures/sky.jpg";
    // std::string bgPath = "/Users/rudolfkischer/Projects/honeycomb/assets/textures/road_sky.hdr";
    // std::string bgPath = "/Users/rudolfkischer/Projects/honeycomb/assets/textures/garden.hdr";
    std::string bgPath = "/Users/rudolfkischer/Projects/honeycomb/assets/textures/rock.hdr";
    // std::string bgPath = "/Users/rudolfkischer/Projects/honeycomb/assets/textures/meadow.hdr";
    GLuint bgTexture = loadTexture(bgPath.c_str(), &bgWidth, &bgHeight, &bgChannels);

    double prevTime = glfwGetTime();

    GLuint rngStateTex;
    glGenTextures(1, &rngStateTex);
    glBindTexture(GL_TEXTURE_2D, rngStateTex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, TEX_WIDTH, TEX_HEIGHT, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, rngStateTex, 0);


    // Render loop
    while (!glfwWindowShouldClose(window))
    {    
        double deltaTime = glfwGetTime() - prevTime;
        prevTime = glfwGetTime();
        GLuint fps = 1 / deltaTime;

        int raySamples = fps / 30;

        // print the ray samples
        // std::cout << "Ray samples: " << raySamples << std::endl;

        updateCameraPosition(window, float(deltaTime));
        glfwGetFramebufferSize(window, &SCR_WIDTH, &SCR_HEIGHT); // update screen size

        bool camMoved = getCamMoved(camPos, prevCamPos, camFront, prevCamFront);
        if (camMoved) { frameCount = 0; } else { frameCount++; }
        fpsCounter(&lastTime, &nbFrames);

        setRenderTarget(shaderProgram, renderTexture->FBO, TEX_WIDTH, TEX_HEIGHT); // render to canvas
        setShaderProgramUniforms(shaderProgram, frameCount, TEX_WIDTH, TEX_HEIGHT, numOfModels, sphereData, camPos, camMoved, modelSize);
        glUniform1fv(glGetUniformLocation(shaderProgram, "materialdata"), numOfMaterials * materialSize, materialData);
        glUniform1i(glGetUniformLocation(shaderProgram, "raySamples"), raySamples);
        //aabb
        glUniform1fv(glGetUniformLocation(shaderProgram, "aabbdata"), numOfAABBs * aabbSize, aabbData);
        // add the number of aabbs
        glUniform1i(glGetUniformLocation(shaderProgram, "numAABBs"), numOfAABBs);

        // bind the rng state texture
        glActiveTexture(GL_TEXTURE4);
        glBindTexture(GL_TEXTURE_2D, rngStateTex);
        glUniform1i(glGetUniformLocation(shaderProgram, "rngState"), 4);


        // pass in the background texture
        glActiveTexture(GL_TEXTURE3);
        glBindTexture(GL_TEXTURE_2D, bgTexture);
        glUniform1i(glGetUniformLocation(shaderProgram, "bgTexture"), 3);
        // pass in bg texture size
        glUniform1i(glGetUniformLocation(shaderProgram, "bgWidth"), bgWidth);
        glUniform1i(glGetUniformLocation(shaderProgram, "bgHeight"), bgHeight);


        glActiveTexture(GL_TEXTURE1); // bind previous frame texture
        glBindTexture(GL_TEXTURE_2D, renderTexture->prevTexture);
        glUniform1i(glGetUniformLocation(shaderProgram, "prevFrame"), 1);
        quad->draw();

        // write the color attachment1 back to the rng state texture
        glBindFramebuffer(GL_READ_FRAMEBUFFER, renderTexture->FBO);
        // glFramebufferTexture2D(GL_READ_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, renderTexture->texture, 0);
        // glBindFramebuffer(GL_DRAW_FRAMEBUFFER, renderTexture->copyFBO);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, rngStateTex, 0);
        // glBlitFramebuffer(0, 0, TEX_WIDTH, TEX_HEIGHT, 0, 0, TEX_WIDTH, TEX_HEIGHT, GL_COLOR_BUFFER_BIT, GL_NEAREST);
        glBindFramebuffer(GL_FRAMEBUFFER, 0);


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
    GLenum drawBuffers[] = {GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1};
    glDrawBuffers(2, drawBuffers);
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

float* crossProduct(float* a, float* b) {
    float* result = new float[3];
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
    return result;
}

void updateCameraPosition(GLFWwindow* window, double deltaTime) {


    // if control is pressed, up the
    if (glfwGetKey(window, GLFW_KEY_LEFT_CONTROL) == GLFW_PRESS)
        speed = 10.0;
    else
        speed = 1.00;
    
    speed *= deltaTime;

    float* up = new float[3]{0.0, 1.0, 0.0};

    float* camRight = crossProduct(camFront, up);


    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) {
        camPos[0] += speed * camFront[0];
        camPos[1] += speed * camFront[1];
        camPos[2] += speed * camFront[2];
    }
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) {
        camPos[0] -= speed * camFront[0];
        camPos[1] -= speed * camFront[1];
        camPos[2] -= speed * camFront[2];
    }

    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) {
        camPos[0] -= speed * camRight[0];
        camPos[1] -= speed * camRight[1];
        camPos[2] -= speed * camRight[2];
    }
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) {
        camPos[0] += speed * camRight[0];
        camPos[1] += speed * camRight[1];
        camPos[2] += speed * camRight[2];
    }
    if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS)
        camPos[1] += speed;
    if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS)
        camPos[1] -= speed;
    
    // turn camera left and right with q and e
    if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS) {
        float angle = -0.01;
        float newCamFront[3];
        newCamFront[0] = camFront[0] * cos(angle) - camFront[2] * sin(angle);
        newCamFront[1] = camFront[1];
        newCamFront[2] = camFront[0] * sin(angle) + camFront[2] * cos(angle);
        for (int i = 0; i < 3; i++) {
            camFront[i] = newCamFront[i];
        }
    }
    if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS) {
        float angle = 0.01;
        float newCamFront[3];
        newCamFront[0] = camFront[0] * cos(angle) - camFront[2] * sin(angle);
        newCamFront[1] = camFront[1];
        newCamFront[2] = camFront[0] * sin(angle) + camFront[2] * cos(angle);
        for (int i = 0; i < 3; i++) {
            camFront[i] = newCamFront[i];
        }
    }
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

bool getCamMoved(float* camPos, float* prevCamPos, float* camFront, float* prevCamFront) {
    bool camMoved = false;
    for (int i = 0; i < 3; i++) {
        if (abs(camPos[i] - prevCamPos[i]) > 0.0001) 
            camMoved = true;
        prevCamPos[i] = camPos[i];
        if (abs(camFront[i] - prevCamFront[i]) > 0.0001)
            camMoved = true;
        prevCamFront[i] = camFront[i];
    }
    return camMoved;
}

void setShaderProgramUniforms(unsigned int shaderProgram, int frameCount, int TEX_WIDTH, int TEX_HEIGHT, int numOfModels, float* sphereData, float* camPos, bool camMoved, int modelSize) {
    glUniform1i(glGetUniformLocation(shaderProgram, "timePassed"), frameCount);
    glUniform1i(glGetUniformLocation(shaderProgram, "texWidth"), TEX_WIDTH);
    glUniform1i(glGetUniformLocation(shaderProgram, "texHeight"), TEX_HEIGHT);
    glUniform1i(glGetUniformLocation(shaderProgram, "numSpheres"), numOfModels);
    glUniform1fv(glGetUniformLocation(shaderProgram, "spheredata"), numOfModels * modelSize, sphereData);
    glUniform3fv(glGetUniformLocation(shaderProgram, "camPos"), 1, camPos);
    glUniform3fv(glGetUniformLocation(shaderProgram, "camFront"), 1, camFront);
    glUniform1i(glGetUniformLocation(shaderProgram, "camMoved"), camMoved);
    glUniform1f(glGetUniformLocation(shaderProgram, "seed1"), dis(gen));
    glUniform1f(glGetUniformLocation(shaderProgram, "uTime"), glfwGetTime());

}



