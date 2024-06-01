#include "Quad.h"



Quad::Quad() {

    //set texture to 0
    texture = 0;

    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);

    // opengl is a state machine and only works with one  vao at a time
    // this will tell opengl to use vao for now
    glBindVertexArray(VAO);
    // this tells opgl to use the specific vbo for now
    // the first argument indicates what kind of buffer
    // it can only bind one , for each kind of buffer
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    // glbufferdata will allocate and optionally initialize the buffer in memory
    // the first agurment, indicates which buffer we are reffering to
    // the second argument indicates the data types size (float64, int32 etc)
    // the third argument is a pointer to the data
    // the last argument indicates how it will be used, allowing for opengl to optimize
    // GL_STATIC_DRAW indicates that it will not be modified over the course of the program
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
    // same does for this buffer array
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

    //This indicates the shape of the individual components of each vertex
    // the first argument indicates the id of the pointer layout, (you can often only have 16-32 of these depending on the hardware)
    // the second argument indicates the length
    // the third is the data type, the fourth tells opengl if itshould normalize the data
    // the fifth argument tells opengl how many bytes it needs to go before it finds the next one
    // the 6th positions indicates the starting offset
    // Position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    // Texture coordinate attribute
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);

    // Unbind VAO by setting bind to 0 (null pointer)
    glBindVertexArray(0);
}

Quad::~Quad() {
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteBuffers(1, &EBO);
}

void Quad::draw() {
    // bind texture
    if (texture) {
      glActiveTexture(GL_TEXTURE0);
      glBindTexture(GL_TEXTURE_2D, texture);
    }

    glBindVertexArray(VAO);

    glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);
}


