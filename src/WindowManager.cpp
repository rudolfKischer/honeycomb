#include "WindowManager.h"
#include <GLFW/glfw3.h>
#include <iostream>



Window::Window(int width, int height, const char* title) {
  // the last two null arguments are used for, indicating:
    // 1. fullscreen (null means windowed) , if you wanted to do fullscreen, you would have to get monitor information, glfwGetMonitors()
    // 2. shared context between windows. If you had already previously made a window, you could pass it in here for it to be used
    this->window = createWindow(width, height, title);
    if (!this->window) { std::cout << "Failed to create window" << std::endl; }
}

Window::~Window() {
    // glfwTerminate();
}


GLFWwindow* Window::createWindow(int width, int height, const char* title) {
    // Initialize the library
    if (!glfwInit())
        return NULL;
    hints();

    // 1. fullscreen (null means windowed) , if you wanted to do fullscreen, you would have to get monitor information, glfwGetMonitors()
    // 2. shared context between windows. If you had already previously made a window, you could pass it in here for it to be used
    GLFWwindow* window = glfwCreateWindow(width, height, title, NULL, NULL);
    
    if (!window) { glfwTerminate(); return NULL; }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(0); // enable vsync
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);


    //ADD KEY CALLBACKS
    // glfwSetKeyCallback(window, key_callback);




    //TODO: ADD MOUSE BUTTON CALLBACKS
    //TODO: ADD CURSOR POSITION CALLBACKS

    return window;
}

void Window::hints() {
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3); // Tell GLFW that it should support OpenGL 3.3 or newer
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE); // mark which functionality it should and shouldnt support
#ifdef __APPLE__ // mark it not to support deprecated functionality because 
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // mac only supports a forward compatible core
#endif
    glfwWindowHint(GLFW_SCALE_TO_MONITOR, GLFW_TRUE);
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height) { 
  glViewport(0, 0, width, height); 
}

// void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
//   if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
//     glfwSetWindowShouldClose(window, true);
//   }
// }

GLFWwindow* make_window(int width, int height, const char* title) {
  Window* window = new Window(width, height, title);
  return window->window;
}



