#include "RenderTexture.h"
// #include <GLFW/glfw3.h>

RenderTexture::RenderTexture(int width, int height) {
    this->width = width;
    this->height = height;
    init();
}

RenderTexture::~RenderTexture() {
    glDeleteFramebuffers(1, &FBO);
    glDeleteTextures(1, &texture);
    glDeleteTextures(1, &prevTexture);
    glDeleteFramebuffers(1, &copyFBO);
}

void RenderTexture::init() {
    glGenFramebuffers(1, &FBO);
    glBindFramebuffer(GL_FRAMEBUFFER, FBO);
    glGenTextures(1, &texture);
    glGenTextures(1, &prevTexture);

    glBindTexture(GL_TEXTURE_2D, texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, texture, 0);

    //init prevTexture
    glBindTexture(GL_TEXTURE_2D, prevTexture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);




    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
        std::cout << "ERROR::FRAMEBUFFER:: Framebuffer is not complete!" << std::endl;
    }
    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    glGenFramebuffers(1, &copyFBO);
}

void RenderTexture::bind() {
    glBindFramebuffer(GL_FRAMEBUFFER, FBO);
    glViewport(0, 0, width, height);
}

void RenderTexture::unbind() {
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void RenderTexture::copy() {
    // Bind the read framebuffer (FBO) and attach the source texture
    glBindFramebuffer(GL_READ_FRAMEBUFFER, FBO);
    glFramebufferTexture2D(GL_READ_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, texture, 0);

    // Bind the draw framebuffer (copyFBO) and attach the destination texture
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, copyFBO);
    glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, prevTexture, 0);

    // Perform the blit (copy operation)
    glBlitFramebuffer(0, 0, width, height, 0, 0, width, height, GL_COLOR_BUFFER_BIT, GL_NEAREST);

    // Unbind the framebuffers
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}


