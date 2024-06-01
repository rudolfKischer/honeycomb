    // bool inShadow(Ray ray, float lightDistance) {
    //     for (int i = 0; i < numSpheres; ++i) {
    //         float t = intersectSphere(ray, spheres[i]);
    //         if (t > 0.0 && t < lightDistance) {
    //             return true; // Shadow ray is blocked
    //         }
    //     }
    //     return false;
    // }

    


    




    //RAY TRACING IMPLEMENTATION

    // float lightDistance = length(lightPos - p);
    // Ray shadowRay;
    // shadowRay.origin = p + normal * 0.0001; // Offset to avoid self-shadowing
    // shadowRay.direction = lightDir;

    // bool shadow = inShadow(shadowRay, lightDistance);

    // // Calculate light contribution
    // float diff = shadow ? 0.0 : max(dot(normal, lightDir), 0.0); // Diffuse light
    // vec3 viewDir = normalize(-ray.direction); // Camera direction
    // vec3 halfwayDir = normalize(lightDir + viewDir);
    // float spec = shadow ? 0.0 : pow(max(dot(normal, halfwayDir), 0.0), 10.0); // Specular highlight

    // vec3 col =  0.6 * closestColor * diff + vec3(1.0) * spec * 0.3; // Combine diffuse and specular
    // col += 0.1 * closestColor; // Add ambient light

    // FragColor = vec4(clamp(col, 0.0, 1.0), 1.0); // Output color with alpha

