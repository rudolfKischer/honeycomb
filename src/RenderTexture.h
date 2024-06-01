#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <iostream>



class RenderTexture {
  public:
      RenderTexture(int width, int height);
      ~RenderTexture();

      void bind();

      void unbind();

      void copy();

  // private:
      void init();
      GLuint FBO, texture, prevTexture, copyFBO;
      GLuint width, height;

};

