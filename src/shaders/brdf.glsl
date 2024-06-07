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


float D_ggx(vec3 h, vec3 n, float alpha) {
    float alpha2 = alpha *alpha;
    float nh = dot(n, h);
    float nh2 = nh * nh;
    float denominator = nh2 * (alpha2 - 1.0) + 1.0;
    return alpha2 * invPI / (denominator * denominator);
}

float G_ggx_shlick(vec3 v, vec3 n, float alpha) {
  float k = alpha * alpha * 0.5;
  float NoV = max(dot(n, v), 0.0);
  return NoV / (NoV * (1.0 - k) + k);

}

float G_smith(vec3 n, vec3 v, vec3 l, float alpha) {
    return G_ggx_shlick(v, n, alpha) * G_ggx_shlick(l, n, alpha);
}

float G_implicit(vec3 n, vec3 v, vec3 l) {
    float NoV = max(dot(n, v), 0.0);
    float NoL = max(dot(n, l), 0.0);
    return NoV * NoL;
}

float G_cook_torrance(vec3 n, vec3 v, vec3 l) {
    vec3 h = normalize(v + l);

    float NoV = dot(n, v);
    float NoL = dot(n, l);
    float NoH = dot(n, h);
    float VoH = dot(v, h);

    float coeff = 2.0 * (NoH) / VoH;
    return min(1.0, min(coeff * NoV, coeff * NoL));
}

vec3 F_shlick(vec3 h, vec3 v, vec3 F0) {
  return F0 + (1.0 - F0) * pow(1.0 - dot(h, v), 5.0);
}

vec3 specular(vec3 h, vec3 n, vec3 v, vec3 l, float alpha, vec3 F0) {
  float D = D_ggx(h, n, alpha);
  float G = G_smith(n, v, l, alpha);
  vec3 F = F_shlick(h, v, F0);
  float NoV = max(dot(n, v), 0.0);
  float NoL = max(dot(n, l), 0.0);
  vec3 specular = D * G * F / (4.0 * NoV * NoL);
  return specular;    
}

vec3 diffuse(vec3 F0, vec3 albedo, vec3 h, vec3 v) {
  vec3 fresnel = F_shlick(h, v, F0);
  return albedo * (1.0 - fresnel);
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

vec3 light(Hit hit, inout vec3 T, inout vec3 L) {
  vec3 n = hit.n;
  vec3 v = hit.ray.direction;
  vec3 h = randomInUnitHemisphere(gSeed, n);
  vec3 l = reflect(-v, h);
  float alpha = hit.mat.orm.y;
  float metal = hit.mat.orm.z;
  vec3 F0_metal = hit.mat.color;
  vec3 F0_dielectric = vec3(0.4);


  vec3 brdf;
  // if the material is metallic, it can either have a metallic specular reflection, or a diffuse reflection, with the probability based on the metallic value
  // if the material is dielectric, it can either have a specular reflection, or a diffuse reflection, with 50% probability
  // in reality, both occur simultaneously, but instead we simulate it with multiple rays
  
  // the interaction, is either a metallic interaction, or a dielectric interaction, with probability based on the metallic value
  bool metallicInteraction = hash1(gSeed) < metal;
  bool isDiffuseInteraction = hash1(gSeed) < 0.5;
  brdf = diffuse(F0_dielectric, hit.mat.color, h, v);
  // if (isDiffuseInteraction) {
  // } else {
    // brdf = specular(h, n, v, l, alpha, F0_dielectric);
  // }
  // if (metallicInteraction) {
  //   brdf = specular(h, n, v, l, alpha, F0_metal);
  // } else {
  //   if (isDiffuseInteraction) {
  //     brdf = diffuse(F0_dielectric, hit.mat.color);
  //   } else {
  //     brdf = specular(h, n, v, l, alpha, F0_dielectric);
  //   }
  // }

  L += T * hit.mat.emission;
  T *= hit.mat.color;

  return l;
}
