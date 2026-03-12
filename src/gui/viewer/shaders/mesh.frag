#version 450

layout(location = 0) in vec3 fragNormal;
layout(location = 1) in vec3 fragWorldPos;
layout(location = 2) in float fragQuality;

layout(location = 0) out vec4 outColor;

layout(set = 0, binding = 0) uniform SceneUBO {
    vec3 lightDir;
    float _pad0;
    vec3 eyePos;
    float _pad1;
    vec4 baseColor;
    int colorMode;  // 0=solid, 1=angle error, 2=convexity, 3=degree
} scene;

// Heat map: blue (good) -> green -> yellow -> red (bad)
vec3 heatMap(float t) {
    t = clamp(t, 0.0, 1.0);
    if (t < 0.25) return mix(vec3(0.0, 0.0, 1.0), vec3(0.0, 1.0, 1.0), t * 4.0);
    if (t < 0.50) return mix(vec3(0.0, 1.0, 1.0), vec3(0.0, 1.0, 0.0), (t - 0.25) * 4.0);
    if (t < 0.75) return mix(vec3(0.0, 1.0, 0.0), vec3(1.0, 1.0, 0.0), (t - 0.50) * 4.0);
    return mix(vec3(1.0, 1.0, 0.0), vec3(1.0, 0.0, 0.0), (t - 0.75) * 4.0);
}

void main() {
    vec3 N = normalize(fragNormal);
    vec3 L = normalize(scene.lightDir);
    vec3 V = normalize(scene.eyePos - fragWorldPos);
    vec3 H = normalize(L + V);

    // Two-sided lighting: flip normal if facing away from camera
    if (dot(N, V) < 0.0) N = -N;

    float ambient = 0.15;
    float diffuse = max(dot(N, L), 0.0) * 0.65;
    float specular = pow(max(dot(N, H), 0.0), 32.0) * 0.3;

    vec3 color;
    if (scene.colorMode == 0) {
        color = scene.baseColor.rgb;
    } else {
        color = heatMap(fragQuality);
    }

    vec3 lit = color * (ambient + diffuse) + vec3(specular);
    outColor = vec4(lit, 1.0);
}
