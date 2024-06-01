#version 330 core

#define MAX_SPHERES 20
const int sphereDatum = 13;
out vec4 FragColor;
in vec2 TexCoord;
uniform sampler2D prevFrame;
uniform sampler2D bgTexture;
uniform int bgWidth;
uniform int bgHeight;
uniform int timePassed;
uniform int texWidth;
uniform int texHeight;
uniform int numSpheres;
uniform bool camMoved;
uniform vec3 camPos;
uniform vec3 camFront;
uniform float spheredata[ MAX_SPHERES * sphereDatum ];
uniform float seed1;
uniform float uTime;
// uniform sampler2D noiseTexture1;

struct Sphere {
    vec3 center;
    float radius;
    vec3 color; // Color attribute of spheres
    vec3 emission; // Emission attribute of spheres
    vec3 orm;
};

struct Ray {
    vec3 origin;
    vec3 direction;
};


const int RAY_BOUNCE_LIMIT =10;
const int RAY_SAMPLES = 1;
float gSeed = 0.0;
Sphere spheres[MAX_SPHERES];


const float epsilon = 0.05;
const float PI = 3.14159265359;
const float invPI = 1.0 / PI;

uint baseHash(uvec2 p)
{
    p = 1103515245U * ((p >> 1U) ^ (p.yx));
    uint h32 = 1103515245U * ((p.x) ^ (p.y>>3U));
    return h32 ^ (h32 >> 16);
}

vec3 hash3(inout float seed)
{
    uint n = baseHash(floatBitsToUint(vec2(seed += 0.1, seed += 0.1)));
    uvec3 rz = uvec3(n, n * 16807U, n * 48271U);
    return vec3(rz & uvec3(0x7fffffffU)) / float(0x7fffffff);
}

float hash1(inout float seed)
{
    uint n = baseHash(floatBitsToUint(vec2(seed += 0.1, seed += 0.1)));
    return float(n & 0x7fffffffU) / float(0x7fffffff);

}

vec3 randomInUnitSphere(inout float seed)
{   
    vec3 h = hash3(seed) * vec3(2.0, 6.28318530718, 1.0) - vec3(1.0, 0.0, 0.0) ;
    float phi = h.y;
    float r = pow(h.z, 1.0/3.0);
	return r * vec3(sqrt(1.0 - h.x * h.x) * vec2(sin(phi), cos(phi)), h.x);
}


// Sphere structure

// lr, lg, lb: emission color
// x, y, z, r, R, G, B, lr, lg, lb
// 0, 1, 2, 3, 4, 5, 6, 7, 8, 9
// sphere unpacker function
Sphere unpackSphere(int index) {
    int offset = index * sphereDatum;
    vec3 center = vec3(spheredata[offset], spheredata[offset + 1], spheredata[offset + 2]);
    float radius = spheredata[offset + 3];
    vec3 color = vec3(spheredata[offset + 4], spheredata[offset + 5], spheredata[offset + 6]);
    vec3 emission = vec3(spheredata[offset + 7], spheredata[offset + 8], spheredata[offset + 9]);
    vec3 orm = vec3(spheredata[offset + 10], spheredata[offset + 11], spheredata[offset + 12]);
    return Sphere(center, radius, color, emission, orm);
}

// unpack the spheres


// Intersect ray with a sphere, return distance, -1 if no hit
float intersectSphere(Ray ray, Sphere sphere) {
    vec3 oc = ray.origin - sphere.center;
    float b = dot(oc, ray.direction);
    float c = dot(oc, oc) - sphere.radius * sphere.radius;
    float h = b * b - c;
    if (h < 0.0) return -1.0; // No intersection
    return -b - sqrt(h);
}

void trace(Ray ray, out float t, out int sphere, int prevSphere) {
    float minT = 99999.0;
    int closestIndex = -1;
    for (int i = 0; i < numSpheres; ++i) {
        if (i == prevSphere) continue;
        float t = intersectSphere(ray, spheres[i]);
        if (t > 0.0 && t < minT) {
            minT = t;
            closestIndex = i;
        }
    }
    t = minT;
    sphere = closestIndex;
}

vec3 f_schlick(float cosTheta, vec3 f0) {
    return f0 + (vec3(1.0) - f0) * pow(1.0 - cosTheta, 5.0);
}

float D_GGX(in float NoH, in float roughness) {
    float a = roughness * roughness;
    float a2 = a * a;
    float NoH2 = NoH * NoH;
    float b = (NoH2 * (a2 - 1.0) + 1.0);
    return a2 * invPI / (b * b);
}

float G1_GGX_smith(in float NoV, in float roughness) {
    float a = roughness * roughness;
    float k = a / 2.0;
    return max(NoV, 0.001) / (NoV * (1.0 - k) + k);
}

float G_Smith(in float NoV, in float NoL, in float roughness) {
    return G1_GGX_smith(NoL, roughness) * G1_GGX_smith(NoV, roughness);
}

float f_shlick90(float cosineTheta, float f0, float f90) {
    return f0 + (f90 - f0) * pow(1.0 - cosineTheta, 5.0);
}

float disneyDiffuseFactor(float NoV, float NoL, float VoH, float roughness) {
    float alpha = roughness * roughness;
    float F90 = 0.5 * VoH * VoH * alpha + 0.25;
    float F_in = f_shlick90(NoL, 1.0, F90);
    float F_out = f_shlick90(NoV, 1.0, F90);
    return F_out * F_in;
}

vec3 brdfmicrofacet(vec3 l, vec3 v, in vec3 n, in vec3 albedo, in float metallic, in float roughness, in float reflectance) {
    roughness = max(roughness, 0.001);

    vec3 h = normalize(v + l);

    float NoV = max(dot(n, v), 0.0);
    float NoL = max(dot(n, l), 0.0);
    float NoH = max(dot(n, h), 0.0);
    float VoH = max(dot(v, h), 0.0);  

    float reflect = max(1-roughness, 0.001);


    vec3 f_0 = vec3(0.16 * reflect * reflect);
    f_0 = mix(f_0, albedo, metallic);
  
    vec3 F = normalize(f_schlick(NoV, f_0));

    float D = D_GGX(NoH, roughness);
    float G = G_Smith(NoL, NoV, roughness);

    vec3 specular = (F * D * G) / (4.0 * max(0.001, NoL) * max(0.001, NoV));

    vec3 rhoD = albedo;

    rhoD *= 1.0 - metallic;

    vec3 diffuse = rhoD * disneyDiffuseFactor(NoV, NoL, VoH, roughness);

    return diffuse + specular;

}

vec2 scaleToFitTexCoords(int width, int height) {
    vec2 screenCoords = gl_FragCoord.xy / vec2(texWidth, texHeight);
    vec2 corrected;

    float textureAspect = float(width) / float(height);
    float screenAspect = float(texWidth) / float(texHeight);

    if (textureAspect > screenAspect) {
        float scale = screenAspect / textureAspect;
        corrected = vec2(screenCoords.x, screenCoords.y * scale + (1.0 - scale) * 0.5);
    } else {
        float scale = textureAspect / screenAspect;
        corrected = vec2(screenCoords.x * scale + (1.0 - scale) * 0.5, screenCoords.y);
    }

    corrected.y = 1.0 - corrected.y;

    return corrected;
}



void main()
{   
    vec2 uv = (gl_FragCoord.xy - 0.5 * vec2(texWidth, texHeight)) / texHeight;
    vec4 prevPixelColor = texture(prevFrame, uv - 0.5);

    // output the background texture as a test
    // vec2 bguv = (gl_FragCoord.xy - 0.5 * vec2(bgWidth, bgHeight)) / bgHeight;
    // bguv.y = -bguv.y;
    // FragColor = texture(bgTexture, scaleToFitTexCoords(bgWidth, bgHeight));
    // return;

    gSeed = float(baseHash(floatBitsToUint(gl_FragCoord.xy))) / float(0xffffffffU) + uTime + seed1 * uTime;

    for (int i = 0; i < numSpheres; ++i) {
        spheres[i] = unpackSphere(i);
    }

    // Camera setup
    float viewportHeight = 2.0;
    float viewportWidth = float(texWidth) / float(texHeight) * viewportHeight;
    float focalLength = 1.0;


    // Create ray from camera to pixel

    vec3 camUp = vec3(0.0, 1.0, 0.0);
    vec3 camRight = normalize(cross(camFront, camUp));

  
    
    // make sure to include camFront direction
    vec3 startRay = normalize(vec3(uv.x * viewportWidth * camRight + uv.y * viewportHeight * camUp + focalLength * camFront));





    vec4 pixelColor = vec4(0.0, 0.0, 0.0, 1.0);

    vec3 T; // throughput, This number represents the current "color" of your filter through which light passes
    vec3 L; // Radiance, this represents the accumulated light, that we add to upon reaching a light
    Ray ray;
    float rayContrib = 1.0 / float(RAY_SAMPLES);
    for (int s=0; s<RAY_SAMPLES; s++) {
        // reset T and L values
        T = vec3(1.0);
        L = vec3(0.0);
        ray.origin = camPos ;
        ray.direction = startRay ;//+ randomInUnitSphere(gSeed) * 0.0005;
        int prevSphere = -1;

        for (int i = 0; i < RAY_BOUNCE_LIMIT; ++i) {
          int closestIndex;
          float minT;

          // Check for intersection with each sphere
          trace(ray, minT, closestIndex, prevSphere);

          if (closestIndex == -1) {
              // vec3 backgroundColor = vec3(1.0, 1.0, 1.0);

              // sample with spherical coordinates
              vec3 dir = normalize(ray.direction);
              float phi = atan(dir.z, dir.x);
              float theta = acos(dir.y);
              vec3 backgroundColor = texture(bgTexture, vec2((phi + PI) / (2.0 * PI), theta / PI)).xyz;
              L += T * backgroundColor;// * max(dot(ray.direction, lightDir), 0.0);
              break;
          }

          // if emission is present, we break
          L += T * spheres[closestIndex].emission;
          if (length(spheres[closestIndex].emission) > 0.5) {
              break;
          }

          //albedo color
          vec3 closestColor = spheres[closestIndex].color; 
          float opacity = spheres[closestIndex].orm.x;

          vec3 p = ray.origin + minT * ray.direction; // Point of intersection
          vec3 normal = normalize(p - spheres[closestIndex].center);


          vec3 reflection = reflect(ray.direction, normal);
          vec3 inDirection = ray.direction;
          vec3 outDirection;

          // create a new ray
          outDirection = randomInUnitSphere(gSeed);

          // outDirection = normalize(outDirection + normal);

          float roughness = spheres[closestIndex].orm.y;
          float metallic = spheres[closestIndex].orm.z;


          // metal reflection
          // outDirection = normalize(reflection * (metallic) * (1.0 - roughness) + randomInUnitSphere(gSeed) * roughness);

          // dielectric reflection using fresnel
          // we use the fresnel coefficient to determine the reflection
          float air = 1.0;
          float glass = 1.4;
          float R_0 = pow((air - glass) / (air + glass), 2);
          float R = R_0 + (1 - R_0) * pow(1 - dot(-ray.direction, normal), 3.0);

          // FragColor = vec4(vec3(R), 1.0);
          // return;

          // FragColor = vec4(vec3(R), 1.0);
          vec3 dialectricReflection = 2 * reflection * R * (1-metallic) * (1.0 - roughness);
          vec3 metalReflection = reflection * (1.0 - R) * (1.0 - roughness) * metallic;


          outDirection = normalize(randomInUnitSphere(gSeed) * (1.0 - R) * roughness + dialectricReflection + metalReflection);

          float transmit = hash1(gSeed);
          if (opacity > 0.0) {
              float ior = 1.4;

              // light should get bent towards the normal

              

              float reflect = hash1(gSeed);
              vec3 randDir = randomInUnitSphere(gSeed);
              if (reflect < R) {
                  if (dot(randDir, normal) > 0.0) {
                      randDir = -randDir;
                  }
                  outDirection = reflection + randDir * roughness;
              } else {
                  if (dot(randDir, normal) < 0.0) {
                      randDir = -randDir;
                  }
                  outDirection = refract(inDirection, normal, 1.0 / ior) + randDir * roughness;
              }
          } 



          ray.direction = outDirection;
          ray.origin = p + ray.direction * 0.0001;
          prevSphere = closestIndex;


          if (transmit > opacity) {
              if (dot(ray.direction, normal) < 0.0) {
                  ray.direction = -ray.direction;
              } 
               //* dot(ray.direction, normal) * 2;
              // T = T * brdfmicrofacet(ray.direction, 
              // -inDirection, 
              // normal, 
              // closestColor, 
              // spheres[closestIndex].orm.z, 
              // spheres[closestIndex].orm.y, 
              // 1.00);
              T = T * closestColor;
          }
          


        }
        pixelColor += vec4(L * rayContrib, 1.0) ;

    }


    


    vec4 finalPixelColor = vec4(clamp(pixelColor.xyz, 0.0, 1.0), 1.0);

    float p = 1.0 / float(timePassed);

    // priortize the new frame, if the camera moved
    if (camMoved) {
        p = 0.2;
    }
    
    finalPixelColor = mix(prevPixelColor, finalPixelColor, p);
    FragColor = finalPixelColor;
}
