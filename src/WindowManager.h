
#include <GLFW/glfw3.h>
#include <iostream>


class Window
{
public:
    Window(int width, int height, const char* title);
    ~Window();
    GLFWwindow* window;

private:

    GLFWwindow* createWindow(int width, int height, const char* title);
    void hints();
};

void framebuffer_size_callback(GLFWwindow* window, int width, int height);

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);

GLFWwindow* make_window(int width, int height, const char* title);