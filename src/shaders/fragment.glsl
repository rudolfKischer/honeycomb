#version 330 core

#define MAX_SPHERES 20
#define MAX_MATERIALS 10
#define MAX_AABB 20
#define MAX_MESH 20
const int sphereDatum = 5;
const int materialDatum = 9;
const int aabbDatum = 7;
const int meshDatum = 14;

layout(location = 0) out vec4 FragColor;
layout(location = 1) out vec4 DataBuffer;
in vec2 TexCoord;
uniform sampler2D prevFrame;
uniform sampler2D bgTexture;
uniform sampler2D rngState;
uniform sampler2D vertexBuffer;
uniform sampler2D indexBuffer;

uniform int VBOwidth;
uniform int EBOwidth;
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
uniform float meshdata[ MAX_MESH * meshDatum ];
uniform float seed1;
uniform float uTime;

uniform int FRACTED_SAMPLES;

vec3 VERTICES[100];
int VERTICES_COUNT = 0;





// uniform sampler2D noiseTexture1;

struct bvhNode {
    vec3 minp;
    vec3 maxp;
    int left;
    int right;
};



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

// struct Mesh {
//     int numVertices;
//     int numFaces;
//     int otherthing;
//     int materialIndex;
//     int vertexStart;
//     ivec3 faces[100];
// };

// struct Mesh {
//     int numVertices;
//     int numFaces;
//     int materialIndex;
//     int vertexStart;
//     int facesStart;
// };

struct Mesh {
    vec3 position;
    int materialIndex;
    int vertexStart;
    int facesStart;
    int numVertices;
    int numFaces;
    vec3 minp;
    vec3 maxp;
};

struct Hit {
    float t;
    vec3 n;
    vec3 p;
    int object;
    int type;
    Material mat;
    Ray ray;
};

const int RAY_BOUNCE_LIMIT = 3;
int RAY_SAMPLES = 3;



float gSeed = 0.0;
Sphere spheres[MAX_SPHERES];
AABB aabbs[MAX_SPHERES];
Material materials[MAX_MATERIALS];


const int numOfMeshes = 1;
Mesh meshes[numOfMeshes];

const float epsilon = 0.05;
const float PI = 3.14159265359;
const float invPI = 1.0 / PI;

vec3 getVertex(int index) {
    float u = mod(float(index), float(VBOwidth));
    float v = floor(float(index) / float(VBOwidth));
    u = u / float(VBOwidth);
    v = v / float(VBOwidth);
    return texture(vertexBuffer, vec2(u, v)).xyz;
}

ivec3 getFace(int index) {
    float u = mod(float(index), float(EBOwidth));
    float v = floor(float(index) / float(EBOwidth));
    u = u / float(EBOwidth);
    v = v / float(EBOwidth);
    return ivec3(texture(indexBuffer, vec2(u, v)).xyz);
}


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

Mesh unpackMesh(int index) {
    int offset = index * meshDatum;
    vec3 position = vec3(meshdata[offset], meshdata[offset + 1], meshdata[offset + 2]);
    int materialIndex = int(meshdata[offset + 3]);
    int vertexStart = int(meshdata[offset + 4]);
    int facesStart = int(meshdata[offset + 5]);
    int numVertices = int(meshdata[offset + 6]);
    int numFaces = int(meshdata[offset + 7]);
    vec3 minp = vec3(meshdata[offset + 8], meshdata[offset + 9], meshdata[offset + 10]);
    vec3 maxp = vec3(meshdata[offset + 11], meshdata[offset + 12], meshdata[offset + 13]);
    return Mesh(position, materialIndex, vertexStart, facesStart, numVertices, numFaces, minp, maxp);
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

// Polygon Intersection
vec3 triangleNormal(vec3 v0, vec3 v1, vec3 v2) {
    return normalize(cross(v1 - v0, v2 - v0));
}

float intersectTriangle(vec3 A, vec3 B, vec3 C, Ray ray, inout vec3 weights) {

  vec3 AB = B - A;
  vec3 AC = C - A;

  vec3 N = normalize(cross(AB, AC));

  // dot(A, N) , gives us how for along the normal we travel  until we reach the plane
  

  // given a point P on the plane, and the point that the normal intersects, P_0
  // we have that the vector from P to P_0
  // T = (P_0 - P)

  // in order for P to be on the plane, it must be the case that T , is perpindicular to the plane

  // for something to be perpindcicular, we must have that the the dot product, with then normal is 0

  // dot(T, N) = 0

  // for us, we want P to be some where on the rays line
  // so our point P  is P = r.o + t * r.d

  // where t is the parameterized distance from theray origin in the rays direction

  // T = (p_0 - (r.o + t * r.d) = (A - (r.o + t * r.d)

  // 0 = N o T 
  // 0 = N o (A - (r.o + t * r.d)
  // 0 = (NoA) - (No(r.o)) - (N o (t*r.d)
  // 0 = (NoA) - (No(r.o)) - (t*(N o r.d))
  // - (NoA) + (No(r.o)) = -(t*(N o r.d))
  // (NoA) - (No(r.o)) = -(t*(N o r.d))
  // [ No(A - r.o) ] / (N o r.d) = t
  float NoD = dot(N, ray.direction);
  if (abs(NoD) < 0.0001) return -1.0;
  float NoRA = dot(N, A - ray.origin);
  float t = NoRA / NoD;

  // now lets get point P
  vec3 P = ray.origin + ray.direction * t;

  // we need to check if P is inside the triangle

  // we can write P as some weighted  combination of the three vertices in the triangle
  // P = u * A + v * B + w * C
  // This can be written as a matrix equation
  // where we have the vertices as rows in a matrix, and the weights as a column vector
  // [ A B C ] * [ u v w ] = P
  // A = [ A B C ]
  // x = [ u v w ]
  // P = A * x
  // x = P * A^-1
  // |x| = 1
  // |P * A^-1| = 1

  // another way to look at the weighting, is to use the area of the sub triangles
  // given a point inside our triangle
  // we can connect this point to our triangles vertices
  // this will give use three sub triangles
  // Notice, that if our point is very close to vertex B, then the area of the sub triangle which does not contain B will be large
  // if P is on B, then our sub triangle will have area equal to the whole triangle, or 100% of the area
  // also not that if P is on the lin connecting A and C, then the area of the sub triangle will be 0

  // it is the case that the area of the sub triangle is proportional to the distance from the vertex

  // so if we get the area of all sub triangles, and divide by the area of the whole triangle, we will get the weights
  // this means we have
  // ABC (whole triangle
  // ABP (C, sub triangle)
  // BCP (A, sub triangle)
  // CAP (B, sub triangle)
  // notice that the corresponding sub triangle is the one that does not contain the vertex

  // the the triangle formed by two vectors
  // spans the parallelogram formed by the two vectors
  // that means the area of the triangle is half the area of the parallelogram
  // the area of the parallelogram is given by the cross product of the two vectors, (which is the definition of the cross product)
  // so the area of the triangle is half the magnitude of the cross product of the two vectors
  vec3 AP = P - A;
  vec3 BP = P - B;
  vec3 CP = P - C;
  // vec3 AB = B - A;
  vec3 BC = C - B;
  vec3 CA = A - C;

  vec3 ABC = cross(AB, AC);
  vec3 ABP = cross(AB, AP);
  if (dot(N, ABP) < 0.0) return -1.0;
  vec3 BCP = cross(BC, BP);
  if (dot(N, BCP) < 0.0) return -1.0;
  vec3 CAP = cross(CA, CP);
  if (dot(N, CAP) < 0.0) return -1.0;

  float areaABC = length(ABC);
  float invAreaABC = 1.0 / areaABC;
  float areaABP = length(ABP);
  float areaBCP = length(BCP);
  float areaCAP = length(CAP);
  weights = vec3(areaBCP, areaCAP, areaABP) * invAreaABC;





  return t;
}






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


Hit intersectAABB(Ray ray, AABB aabb) {
    vec3 minp = aabb.minp;
    vec3 maxp = aabb.maxp;
  // normalize the ray direction
    ray.direction = normalize(ray.direction);
    vec3 invDir = 1.0 / (ray.direction + 0.0000000001);
    vec3 LB = (minp - ray.origin) * invDir;
    vec3 UB = (maxp - ray.origin) * invDir;
    vec3 AABBmin = min(LB, UB);
    vec3 AABBmax = max(LB, UB);
    float tmin = minC(AABBmax);
    float tmax = maxC(AABBmin);
    float t;
    if (tmax < 0.0 || tmin < tmax) {
      t = -1.0;
    } else {
      t = tmax;
    }


    Hit hit;
    hit.t = t;
    hit.p = ray.origin + ray.direction * hit.t;
    hit.n = aabbNormal(aabb, hit.p);
    hit.mat = materials[aabb.materialIndex];
    hit.object = 0;
    hit.type = 1;
    hit.ray = ray;
    return hit;
}


// Intersect ray with a sphere, return distance, -1 if no hit
Hit intersectSphere(Ray ray, Sphere sphere) {
    vec3 oc = ray.origin - sphere.center;
    float b = dot(oc, ray.direction);
    float c = dot(oc, oc) - sphere.radius * sphere.radius;
    float h = b * b - c;
    float t;
    if (h < 0.0) {
      t = -1.0;
    } else {
      t = -b - sqrt(h);
    }


    Hit hit;
    hit.t = t;
    hit.p = ray.origin + ray.direction * hit.t;
    hit.n = normalize(hit.p - sphere.center);
    hit.mat = materials[sphere.materialIndex];
    hit.object = 0;
    hit.type = 0;
    hit.ray = ray;
    return hit;
}




Hit intersectMesh(Ray ray, Mesh mesh) {

   

    Hit hit;
    hit.t = 99999999.0;
    hit.object = -1;

    // before checking the faces, we need to check if the ray intersects the aabb
    AABB hitbox = AABB(mesh.minp + mesh.position, mesh.maxp + mesh.position, mesh.materialIndex);
    Hit aabbHit = intersectAABB(ray, hitbox);
    if (aabbHit.t < 0.0) {
        return hit;
    }


    for (int j = 0; j < mesh.numFaces; ++j) {
        // ivec3 face = mesh.faces[j];
        ivec3 face = getFace(mesh.facesStart + j);
        // vec3 A = VERTICES[face.x];
        // vec3 B = VERTICES[face.y];
        // vec3 C = VERTICES[face.z];
        vec3 A = getVertex(mesh.vertexStart + face.x);
        vec3 B = getVertex(mesh.vertexStart + face.y);
        vec3 C = getVertex(mesh.vertexStart + face.z);
        // add position
        A += mesh.position;
        B += mesh.position;
        C += mesh.position;
        vec3 weights = vec3(0.0, 0.0, 0.0);
        float t = intersectTriangle(A, B, C, ray, weights);
        if (t > 0.00001 && t < hit.t) {
            hit.t = t;
            hit.p = ray.origin + ray.direction * hit.t;
            hit.n = triangleNormal(A, B, C);
            hit.mat = materials[mesh.materialIndex];
            // if the normal is facing the same direction as the ray, then we flip it
            if (dot(hit.n, ray.direction) > 0.0) {
                hit.n = -hit.n;
            }

            // change the mat diffuse color based on the weights
            // we want to make it so that the edges are black
            // to do this using the weights, what we do, is if one weights is less than 0.1, then we make the color black
            float threshold = 0.01;
            hit.mat.color = vec3(1.0, 1.0, 1.0);
            if (weights.x < threshold) {
                hit.mat.color = vec3(0.0, 0.0, 1.0);
                hit.mat.emission = vec3(0.0, 0.0, 1.0);
            }
            if (weights.y < threshold) {
                hit.mat.color = vec3(1.0, 0.0, 0.0);
                hit.mat.emission = vec3(1.0, 0.0, 0.0);
            }
            if (weights.z < threshold) {
                hit.mat.color = vec3(0.0, 1.0, 0.0);
                hit.mat.emission = vec3(0.0, 1.0, 0.0);
            }
            hit.object = 0;
            hit.type = 2;
        }
    }

    if (hit.object == -1) {
      hit.t = -1.0;
    }
    return hit;
}

Hit getIntersection(Ray ray, int prevSphere) {
    // intersect spheres
    Hit hit;
    hit.t = 999999.0;
    hit.object = -1;
    for (int i = 0; i < numSpheres; ++i) {
        Hit newHit = intersectSphere(ray, spheres[i]);
        if (newHit.t > 0.0 && newHit.t < hit.t) {
            hit = newHit;
            hit.object = i;
        }
    }

    // intersect aabbs
    for (int i = 0; i < numAABBs; ++i) {
        Hit newHit = intersectAABB(ray, aabbs[i]);
        if (newHit.t > 0.0 && newHit.t < hit.t) {
            hit = newHit;
            hit.object = i;
        }
    }
  
    for (int i = 0; i < numOfMeshes; ++i) {
        Hit newHit = intersectMesh(ray, meshes[i]);
        if (newHit.t > 0.0 && newHit.t < hit.t) {
            hit = newHit;
            hit.object = i;
        }
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

    // vec3 diffuse = normalize(randomInUnitHemisphere(gSeed, hit.n));
    // T = T * hit.mat.color;

    // return diffuse;

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
        // hsv.z *= 1.3;
        // hsv.z = pow(hsv.z, 3.5);
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
    ray.origin = hit.p + ray.direction * 0.1;
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
    startDir = normalize(startDir + randomCosineInUnitSphere(gSeed) *0.001);
    return Ray(camPos, startDir);
}

vec3 samplePixel(Ray ray) {
    int prevSphere = -1;
    vec3 T = vec3(1.0, 1.0, 1.0);
    vec3 L = vec3(0.0, 0.0, 0.0);
    vec3 color = vec3(0.0, 0.0, 0.0);
    for (int i = 0; i < RAY_BOUNCE_LIMIT; ++i) {
      Hit hit = getIntersection(ray, prevSphere);
      if (hit.object == -1) { missedHit(ray, T, L, i); break; }
      
      L += T * hit.mat.emission;
      if (length(hit.mat.emission) > 0.5) { break; } // if hit light source, return light source color
      accumulate(T, L, ray, hit, prevSphere, i);
      color = hit.mat.color;
    }
    return L;
}

vec3 getPixelColor(vec2 uv) {
  vec3 pixelColor = vec3(0.0);
  float sampleContribution = 1.0 / float(RAY_SAMPLES);
  // uv = fract(uv * 2.0);
  // Ray ray = startRay(fract(uv * 2.0));
  // vec2 uvNormalized = (uv + 1) ;
  // uvNormalized = fract(uvNormalized);
  // vec2 fractedShiftedUv = uvNormalized - 1.0;
  vec2 fragpos = (gl_FragCoord.xy * FRACTED_SAMPLES);
  fragpos.x = int(fragpos.x) % texHeight;
  fragpos.y = int(fragpos.y) % texHeight;
  vec2 uv_frac = ((fragpos) - (0.5 * vec2(texWidth, texHeight))) / texHeight;
  Ray ray = startRay(uv_frac);
  // Ray ray = startRay(uv);
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

    vec2 uv = (gl_FragCoord.xy - (0.5 * vec2(texWidth, texHeight))) / texHeight;

    // Mesh testMesh = getTestTriangle();
    // Mesh testMesh = getBoxMesh();
    // Mesh testMesh = getTetraHedronMesh();
    // Mesh testMesh = getOctaHedron();
    // Mesh testMesh = getDodecahedronMesh();

    // Mesh testMesh;
    // //triangle test mesh
    // testMesh.numVertices = 3;
    // testMesh.numFaces = 1;
    // testMesh.materialIndex = 0;
    // testMesh.vertexStart = 0;
    // testMesh.facesStart = 0;


    // meshes[0] = testMesh;
    
    // we want to render multiple passes at once
    // to do this , we will render the same scene multiple times, on the same frame
    // this will look like having for example, 4 copies of the same scene on the same frame
    // to do this we will modify the uv coordinates

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

    for (int i = 0; i < numOfMeshes; ++i) {
        meshes[i] = unpackMesh(i);
    }



    vec3 pixelColor = getPixelColor(uv);
    DataBuffer = vec4(vec3(gSeed), 1.0);
    FragColor = vec4(pixelColor, 1.0);

    // test accessing the first vertex of the vertex buffer
    // use the vertex as the color
    // vec3 color = getVertex(1);
    // FragColor = vec4(color, 1.0);

    // to test the faces, we will use the first face as the color
    // ivec3 face = getFace(0);
    // vec3 color = vec3(face);
    // FragColor = vec4(color, 1.0);

    // for debug purposes display the vertices in the vertex buffer to the screen
    // we want to use the frags uv position to index into the vertex buffer
    // we want to scale to fit the screen based on the number of vertices
    // vec3 color = getVertex(int(gl_FragCoord.x / 20.0));
    // FragColor = vec4(color, 1.0);

    // for debug purposes display the faces in the face buffer to the screen
    // we want to use the frags uv position to index into the face buffer
    // we want to scale to fit the screen based on the number of faces
    // there are 12 faces, we we want to get our uv coord between 0 and 12
    // float u = gl_FragCoord.x / texWidth;
    // float v = gl_FragCoord.y / texHeight;
    // ivec3 face = getFace(int(u * 12.0));
    // vec3 index = texture(indexBuffer, vec2(u, v)).xyz;

    float face_num = (gl_FragCoord.x / texWidth) * 36.0;
    float u = mod(face_num, EBOwidth);
    float v = floor(face_num / EBOwidth);
    // vec3 face = texture(indexBuffer, vec2(u / EBOwidth, v / EBOwidth)).xyz;
    vec3 face = getFace(int(face_num));







    // convet to vec3
    // vec3 color = vec3(u / 12.0, v / 12.0, 0.0);
    vec3 color = vec3(face);
    // vec3 color = vec3(u);
    // divide by the number of vertices (8)
    color /= 107.0;
    // FragColor = vec4(color, 1.0);

}
