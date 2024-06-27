#version 330 core
out vec4 FragColor;
in vec2 TexCoords;
uniform sampler2D screenTexture;
uniform int FRACTED_SAMPLES;


// const int FRACTED_SAMPLES = 2;

// we take multiple samples at the same time
// to do this we draw the scene  FRACTED_SAMPLES * FRACTED_SAMPLES times
// if FRACTED_SAMPLES = 1, then we can just sample the texCoords
// other wise we need to do an average of all the samples
// To do this we get the corresponding position in the top left corner of the screen
// then we iterate by row and col, for the number of fractedsamples
// to get the position of the sample we add the row and col to the top left corner
// we take our texCoords


void main()
{   
    // FragColor = texture(screenTexture, TexCoords);
    // return;
    // vec4 pixelColor = vec4(0.0);
    // for (int row = 0; row < FRACTED_SAMPLES; ++row)
    // {
    //     for (int col = 0; col < FRACTED_SAMPLES; ++col)
    //     {  
    //         vec2 startPos = TexCoords / float(FRACTED_SAMPLES);
    //         vec2 relativePos = vec2(float(row), float(col)) / float(FRACTED_SAMPLES);
    //         vec2 samplePos = startPos + relativePos;
    //         pixelColor += texture(screenTexture, samplePos);
    //     }
    // }
    // pixelColor /= float(FRACTED_SAMPLES * FRACTED_SAMPLES);

// 
    // FragColor = pixelColor;


    // return;







    // apply FXAA (Fast Approximate Anti-Aliasing)
    vec2 texelSize = 1.0 / vec2(textureSize(screenTexture, 0));
    vec3 rgbNW = texture(screenTexture, TexCoords + (vec2(-1.0, -1.0) * texelSize)).xyz;
    vec3 rgbNE = texture(screenTexture, TexCoords + (vec2(1.0, -1.0) * texelSize)).xyz;
    vec3 rgbSW = texture(screenTexture, TexCoords + (vec2(-1.0, 1.0) * texelSize)).xyz;
    vec3 rgbSE = texture(screenTexture, TexCoords + (vec2(1.0, 1.0) * texelSize)).xyz;
    vec3 rgbM  = texture(screenTexture, TexCoords).xyz;

    vec3 luma = vec3(0.299, 0.587, 0.114);
    float lumaNW = dot(rgbNW, luma);
    float lumaNE = dot(rgbNE, luma);
    float lumaSW = dot(rgbSW, luma);
    float lumaSE = dot(rgbSE, luma);
    float lumaM  = dot(rgbM,  luma);

    float lumaMin = min(lumaM, min(min(min(lumaNW, lumaNE), lumaSW), lumaSE));
    float lumaMax = max(lumaM, max(max(max(lumaNW, lumaNE), lumaSW), lumaSE));

    vec2 dir;
    dir.x = -((lumaNW + lumaNE) - (lumaSW + lumaSE));
    dir.y =  ((lumaNW + lumaSW) - (lumaNE + lumaSE));

    float dirReduce = max((lumaNW + lumaNE + lumaSW + lumaSE) * (0.25 * 0.25), 0.001);
    float rcpDirMin = 1.0 / (min(abs(dir.x), abs(dir.y)) + dirReduce);

    dir = min(vec2(8.0, 8.0), max(vec2(-8.0, -8.0), dir * rcpDirMin)) * texelSize;

    vec3 rgbA = 0.5 * (
        texture(screenTexture, TexCoords + dir * (1.0 / 3.0 - 0.5)).xyz +
        texture(screenTexture, TexCoords + dir * (2.0 / 3.0 - 0.5)).xyz);

    vec3 rgbB = rgbA * 0.5 + 0.25 * (
        texture(screenTexture, TexCoords + dir * -0.5).xyz +
        texture(screenTexture, TexCoords + dir * 0.5).xyz);

    float lumaB = dot(rgbB, luma);

    if (lumaB < lumaMin || lumaB > lumaMax)
        FragColor = vec4(rgbA, 1.0);
    else
        FragColor = vec4(rgbB, 1.0);



}