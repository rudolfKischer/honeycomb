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
