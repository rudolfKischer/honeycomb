
#include <glad/glad.h>
#include <GLFW/glfw3.h>

class Quad {
  public:
      unsigned int VAO, VBO, EBO;
      unsigned int texture;
      Quad();
      ~Quad();
      void draw();
  private:
      void init();
      // verts and indices
      float vertices[20] = {
          // positions         // texture coords
          1.0f,  1.0f, 0.0f,  1.0f, 1.0f,   // top right
          1.0f, -1.0f, 0.0f,  1.0f, 0.0f,   // bottom right
          -1.0f, -1.0f, 0.0f,  0.0f, 0.0f,   // bottom left
          -1.0f,  1.0f, 0.0f,  0.0f, 1.0f    // top left 
      };

      unsigned int indices[6] = {
          0, 1, 3, // first triangle
          1, 2, 3  // second triangle
      };






};