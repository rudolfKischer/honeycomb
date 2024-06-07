#version 330 core

#define MAX_SPHERES 20
#define MAX_MATERIALS 10
#define MAX_AABB 20
const int sphereDatum = 5;
const int materialDatum = 9;
const int aabbDatum = 7;
layout(location = 0) out vec4 FragColor;
layout(location = 1) out vec4 DataBuffer;
in vec2 TexCoord;
uniform sampler2D prevFrame;
uniform sampler2D bgTexture;
uniform sampler2D rngState;
uniform int bgWidth;
uniform int bgHeight;
uniform int timePassed;
uniform int texWidth;
uniform int texHeight;
uniform int numSpheres;
uniform int numAABBs;

uniform int raySamples;
uniform bool camMoved;
uniform vec3 camPos;
uniform vec3 camFront;
uniform float spheredata[ MAX_SPHERES * sphereDatum ];
uniform float materialdata[ MAX_MATERIALS * materialDatum ];
uniform float aabbdata[ MAX_SPHERES * aabbDatum ];
uniform float seed1;
uniform float uTime;

// uniform sampler2D noiseTexture1;

struct Sphere {
    vec3 center;
    float radius;
    int materialIndex;
};

struct AABB {
    vec3 minp;
    vec3 maxp;
    int materialIndex;
};

struct Ray {
    vec3 origin;
    vec3 direction;
};

struct Material {
    vec3 color;
    vec3 emission;
    vec3 orm;
};

struct Mesh {
    int numVertices;
    int numIndices;
    int materialIndex;
    vec3 vertices[100];
    vec3 normals[100];
    vec2 uvs[100];
    int indices[100];
};

struct Hit {
    float t;
    vec3 n;
    vec3 p;
    int object;
    Material mat;
    Ray ray;
};

const int RAY_BOUNCE_LIMIT = 4;
int RAY_SAMPLES = 6;



float gSeed = 0.0;
Sphere spheres[MAX_SPHERES];
AABB aabbs[MAX_SPHERES];
Material materials[MAX_MATERIALS];

const float epsilon = 0.05;
const float PI = 3.14159265359;
const float invPI = 1.0 / PI;


uint baseHash(uvec2 p)
{
    p = 1103515245U * ((p >> 1U) ^ (p.yx));
    uint h32 = 1103515245U * ((p.x) ^ (p.y>>3U));
    return h32 ^ (h32 >> 16);
}

float hash1(inout float seed)
{
    uint n = baseHash(floatBitsToUint(vec2(seed , seed)));

    float f = float(n & 0x7fffffffU) / float(0x7fffffff);
    seed = f;
    return f;

}
vec3 hash3(inout float seed)
{
    // uint n1 = baseHash(floatBitsToUint(vec2(seed, seed)));
    float f1 = hash1(seed);
    float f2 = hash1(seed);
    float f3 = hash1(seed);
    return vec3(f1, f2, f3);
}

float gaussian(inout float gseed) {
    float r1 = hash1(gseed);
    float r2 = hash1(gseed);
    float theta = 2.0 * PI * r1;
    float r = sqrt(-2.0 * log(r2));
    return r * cos(theta);
}


vec3 randomInUnitSphere(inout float seed)
{   
    float g1 = gaussian(seed);
    float g2 = gaussian(seed);
    float g3 = gaussian(seed);
    return normalize(vec3(g1, g2, g3));
}

vec3 randomInUnitHemisphere(inout float seed, vec3 normal) {
    vec3 inUnitSphere = randomInUnitSphere(seed);
    return inUnitSphere * sign(dot(inUnitSphere, normal));
}

vec3 randomCosineInUnitSphere(inout float seed)
{
    float r1 = hash1(seed);
    float r2 = hash1(seed);
    float azimuth = 2.0 * PI * r1;
    float zenith = acos(sqrt(r2));
    return vec3(sin(zenith) * cos(azimuth), sin(zenith) * sin(azimuth), cos(zenith));
}

Sphere unpackSphere(int index) {
    int offset = index * sphereDatum;
    vec3 center = vec3(spheredata[offset], spheredata[offset + 1], spheredata[offset + 2]);
    float radius = spheredata[offset + 3];
    int materialIndex = int(spheredata[offset + 4]);
    return Sphere(center, radius, materialIndex);
}



// Material structure
Material unpackMaterial(int index) {
    int offset = index * 9;
    vec3 color = vec3(materialdata[offset], materialdata[offset + 1], materialdata[offset + 2]);
    vec3 emission = vec3(materialdata[offset + 3], materialdata[offset + 4], materialdata[offset + 5]);
    vec3 orm = vec3(materialdata[offset + 6], materialdata[offset + 7], materialdata[offset + 8]);
    return Material(color, emission, orm);
}

// AABB structure
AABB unpackAABB(int index) {
    int offset = index * 7;
    vec3 minp = vec3(aabbdata[offset], aabbdata[offset + 1], aabbdata[offset + 2]);
    vec3 maxp = vec3(aabbdata[offset + 3], aabbdata[offset + 4], aabbdata[offset + 5]);
    int materialIndex = int(aabbdata[offset + 6]);
    return AABB(minp, maxp, materialIndex);
}

// // Polygon Intersection
// float intersectTriangle(vec3 A, vec3 B, vec3 C, Ray ray) {
//     vec3 AB = B - A;
//     vec3 AC = C - A;
//     vec3 N = cross(AB, AC);


// }


// AABB intersection
float minC(vec3 v) {
    return min(min(v.x, v.y), v.z);
}

float maxC(vec3 v) {
    return max(max(v.x, v.y), v.z);
}

vec3 aabbNormal(AABB aabb, vec3 intersectionPoint) {
    // find the closest face  
    vec3 minp = aabb.minp;
    vec3 maxp = aabb.maxp;

    vec3 diffmax = abs(maxp - intersectionPoint);
    vec3 diffmin = abs(minp - intersectionPoint);

    //float max
    float minValue = float(0xffffffffU);
    float minIndex = -1;
    float minOrMax = 0.0;
    for (int i = 0; i < 3; ++i) {
        if (diffmax[i] < minValue) {
            minValue = diffmax[i];
            minIndex = i;
            minOrMax = 1.0;
        }
        if (diffmin[i] < minValue) {
            minValue = diffmin[i];
            minIndex = i;
            minOrMax = -1.0;
        }
    }
    vec3 normal = vec3(0.0, 0.0, 0.0);
    normal[int(minIndex)] = minOrMax;

    return normal;
}


float intersectAABB(vec3 minp, vec3 maxp, Ray ray) {
  // normalize the ray direction
    ray.direction = normalize(ray.direction);
    vec3 invDir = 1.0 / (ray.direction + 0.0000000001);
    vec3 LB = (minp - ray.origin) * invDir;
    vec3 UB = (maxp - ray.origin) * invDir;
    vec3 AABBmin = min(LB, UB);
    vec3 AABBmax = max(LB, UB);
    float tmin = minC(AABBmax);
    float tmax = maxC(AABBmin);
    if (tmax < 0.0 || tmin < tmax) return -1.0;
    return tmax;
}


// Intersect ray with a sphere, return distance, -1 if no hit
float intersectSphere(Ray ray, Sphere sphere) {
    vec3 oc = ray.origin - sphere.center;
    float b = dot(oc, ray.direction);
    float c = dot(oc, oc) - sphere.radius * sphere.radius;
    float h = b * b - c;
    if (h < 0.0) return -1.0; // No intersection
    return -b - sqrt(h);
}

Hit getIntersection(Ray ray, int prevSphere) {
    float minT = 99999.0;
    // intersect spheres
    Hit hit;
    hit.object = -1;
    int hitType = -1;
    for (int i = 0; i < numSpheres; ++i) {
        float t = intersectSphere(ray, spheres[i]);
        if (t > 0.0 && t < minT) {
            minT = t;
            hit.object = i;
            hitType = 0;
        }
    }
    // intersect aabbs
    for (int i = 0; i < numAABBs; ++i) {
        float t = intersectAABB(aabbs[i].minp, aabbs[i].maxp, ray);
        if (t > 0.0 && t < minT) {
            minT = t;
            hit.object = i;
            hitType = 1;
        }
    }

    if (hit.object == -1) return hit;
    hit.t = minT;
    hit.p = ray.origin + ray.direction * hit.t;
    if (hitType == 0) {
        hit.n = normalize(hit.p - spheres[hit.object].center);
        hit.mat = materials[spheres[hit.object].materialIndex];
    } else {
        hit.n = aabbNormal(aabbs[hit.object], hit.p);
        hit.mat = materials[aabbs[hit.object].materialIndex];
    }
    return hit;
}

float shlickFresnel(float n1, float n2, vec3 normal, vec3 incident) {
    float r0 = pow((n1 - n2) / (n1 + n2), 2);
    float cosTheta = clamp(dot(-incident, normal), 0.0, 1.0);
    float r = r0 + (1 - r0) * pow(1.0 - cosTheta, 5);
    return r;
}

vec3 D_GGX_sample(float alpha, vec3 n) {
  float theta = atan(sqrt(-alpha * alpha * log(1.0 - hash1(gSeed))));
  float phi = 2.0 * PI * hash1(gSeed);
  vec3 H = vec3(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));

  // transform H to be in the same hemisphere as n
  if (dot(H, n) < 0.0) {
    H = -H;
  }
  return H;
}

float G_ggx_shlick(vec3 v, vec3 n, float alpha) {
  float k = alpha * alpha * 0.5;
  float NoV = max(dot(n, v), 0.0);
  return NoV / (NoV * (1.0 - k) + k);

}

float G_smith(vec3 n, vec3 v, vec3 l, float alpha) {
    return G_ggx_shlick(v, n, alpha) * G_ggx_shlick(l, n, alpha);
}



vec3 sampleOutRay(Ray r_i, Hit hit, inout vec3 T) {

    float diceRoll = hash1(gSeed);
    float diceRoll2 = hash1(gSeed);// I DONT KNOW WHY BUT THESE NEED TO HAPPEN BEFORE THE BRANCHING, otherwise things dont stay random

    float epsilon = 0.01;

    // Sample the microfacet normal
    float surfaceRoughness = hit.mat.orm.y; 
    // vec3 randNormHemisphere = randomInUnitHemisphere(gSeed, hit.n);
    vec3 randNormHemisphere = D_GGX_sample(surfaceRoughness, hit.n);
    vec3 micro_normal = mix(hit.n, randNormHemisphere, surfaceRoughness);


    // shadowing and masking
    if (-dot(micro_normal, r_i.direction) < 0.0) {
      vec3 l = -reflect(r_i.direction, micro_normal);
      float geometry = G_smith(hit.n, -r_i.direction, l, surfaceRoughness);
      T = T * sqrt(geometry);
    }

    // Metallic case
    vec3 reflection = normalize(reflect(r_i.direction, micro_normal));
    float metallic = hit.mat.orm.z;
    if (metallic > epsilon) { 
      T = T * hit.mat.color;
      return reflection; 
    }

    // Dielectric case
    float ior = 1.4;
    float fresnel = shlickFresnel(1.0, ior, hit.n, r_i.direction);
    if (diceRoll < fresnel * (1 - surfaceRoughness)) {
      return reflection;
    }

    // Refracted case
    float opacity = hit.mat.orm.x; 
    vec3 refracted = refract(r_i.direction, micro_normal, 1.0 / ior);
    if (diceRoll2 < opacity) { 
      T = T * hit.mat.color;
      return refracted; 
    }

    // diffuse case
    vec3 diffuse = normalize(randomInUnitHemisphere(gSeed, hit.n));
    T = T * hit.mat.color;

    return diffuse;

}

vec3 rgbToHSV(vec3 color) {
    float cmax = max(color.r, max(color.g, color.b));
    float cmin = min(color.r, min(color.g, color.b));
    float delta = cmax - cmin;
    float hue = 0.0;
    if (delta != 0.0) {
        if (cmax == color.r) {
            hue = (color.g - color.b) / delta;
        } else if (cmax == color.g) {
            hue = 2.0 + (color.b - color.r) / delta;
        } else {
            hue = 4.0 + (color.r - color.g) / delta;
        }
    }
    hue *= 60.0;
    if (hue < 0.0) {
        hue += 360.0;
    }
    float saturation = cmax == 0.0 ? 0.0 : delta / cmax;
    float value = cmax;
    return vec3(hue, saturation, value);
}

vec3 hsvToRGB(vec3 color) {
    float c = color.z * color.y;
    float x = c * (1.0 - abs(mod(color.x / 60.0, 2.0) - 1.0));
    float m = color.z - c;
    vec3 primeColor;
    if (color.x >= 0.0 && color.x < 60.0) {
        primeColor = vec3(c, x, 0.0);
    } else if (color.x >= 60.0 && color.x < 120.0) {
        primeColor = vec3(x, c, 0.0);
    } else if (color.x >= 120.0 && color.x < 180.0) {
        primeColor = vec3(0.0, c, x);
    } else if (color.x >= 180.0 && color.x < 240.0) {
        primeColor = vec3(0.0, x, c);
    } else if (color.x >= 240.0 && color.x < 300.0) {
        primeColor = vec3(x, 0.0, c);
    } else {
        primeColor = vec3(c, 0.0, x);
    }
    return primeColor + vec3(m);
}

void missedHit(Ray ray, inout vec3 T, inout vec3 L, int depth) {
    vec3 dir = normalize(ray.direction);
    float phi = atan(dir.z, dir.x);
    float theta = acos(dir.y);
    vec3 backgroundColor = texture(bgTexture, vec2((phi + PI) / (2.0 * PI), theta / PI)).xyz;
    backgroundColor = sqrt(backgroundColor * 0.25);
    if (depth != 0) {
        // coonvert to hsv, and bring the luminance up by taking it to a power
        // then convert back to rgb
        vec3 hsv = rgbToHSV(backgroundColor);
        hsv.z *= 1.3;
        hsv.z = pow(hsv.z, 1.5);
        backgroundColor = hsvToRGB(hsv);

        // divide by the largest value in the texture to normalize
        // we need to check the whole texture to find the largest value
        L += T * backgroundColor;


        return;
    }
    L += T * backgroundColor;
    return;
}

void accumulate(inout vec3 T, inout vec3 L, inout Ray ray, Hit hit, inout int prevSphere, in int depth) {
    
    // if (length(hit.mat.emission) > 0.5) { return; }
    ray.direction = sampleOutRay(ray, hit, T);
    ray.origin = hit.p + ray.direction * 0.00001;
    prevSphere = hit.object;

    return;
}


void genSeed(inout float genSeed) {
    vec2 uv = (gl_FragCoord.xy - 0.5 * vec2(texWidth, texHeight)) / texHeight;
    vec4 rng = texture(rngState, uv);
    uvec2 fragSeed = uvec2(floatBitsToUint(gl_FragCoord.x), floatBitsToUint(gl_FragCoord.y));
    uint tseed = floatBitsToUint(float(baseHash(uvec2(uint(timePassed), uint(timePassed)))));
    tseed ^= floatBitsToUint(seed1);
    fragSeed.x ^= tseed;
    fragSeed.y ^= floatBitsToUint(float(baseHash(uvec2(tseed,tseed))));
    fragSeed.x ^= floatBitsToUint(seed1);
    fragSeed.y ^= floatBitsToUint(seed1);
    genSeed = float(baseHash(fragSeed)) / float(0xffffffffU);
}

Ray startRay(in vec2 uv) {
    float viewportHeight = 2.0;
    float viewportWidth = float(texWidth) / float(texHeight) * viewportHeight;
    float focalLength = 1.0;
    vec3 camUp = vec3(0.0, 1.0, 0.0);
    vec3 camRight = normalize(cross(camFront, camUp));
    vec3 startDir = normalize(vec3(uv.x * viewportWidth * camRight + uv.y * viewportHeight * camUp + focalLength * camFront));
    startDir = normalize(startDir + randomCosineInUnitSphere(gSeed) *0.00001);
    return Ray(camPos, startDir);
}

vec3 samplePixel(Ray ray) {
    int prevSphere = -1;
    vec3 T = vec3(1.0, 1.0, 1.0);
    vec3 L = vec3(0.0, 0.0, 0.0);
    for (int i = 0; i < RAY_BOUNCE_LIMIT; ++i) {
      Hit hit = getIntersection(ray, prevSphere);
      // return (vec3(hit.n) + 1.0) * 0.5;
      if (hit.object == -1) { missedHit(ray, T, L, i); break; }
      L += T * hit.mat.emission;
      if (length(hit.mat.emission) > 0.5) { break; } // if hit light source, return light source color
      accumulate(T, L, ray, hit, prevSphere, i);
    }
    return L;
}

vec3 getPixelColor(vec2 uv) {
  vec3 pixelColor = vec3(0.0);
  float sampleContribution = 1.0 / float(RAY_SAMPLES);
  Ray ray = startRay(uv);
  for (int s=0; s<RAY_SAMPLES; s++) {
      pixelColor += samplePixel(ray);
  }
  pixelColor *= sampleContribution;
  pixelColor = clamp(pixelColor.xyz, 0.0, 1.0);

  float p = 1.0 / float(timePassed);
  if (camMoved) {
      p = 0.3;
  }
  vec3 prevPixelColor = texture(prevFrame, uv - 0.5).xyz;
  pixelColor = mix(prevPixelColor, pixelColor, p);
  return pixelColor;
}

void main()
{   
    vec2 uv = (gl_FragCoord.xy - 0.5 * vec2(texWidth, texHeight)) / texHeight;
    
    genSeed(gSeed);

    for (int i = 0; i < numSpheres; ++i) {
        spheres[i] = unpackSphere(i);
    }

    for (int i = 0; i < MAX_MATERIALS; ++i) {
        materials[i] = unpackMaterial(i);
    }

    for (int i = 0; i < numAABBs; ++i) {
        aabbs[i] = unpackAABB(i);
    }



    vec3 pixelColor = getPixelColor(uv);
    DataBuffer = vec4(vec3(gSeed), 1.0);
    FragColor = vec4(pixelColor, 1.0);
}
