#version 450

layout(push_constant) uniform PushConstants {
    mat4 mvp;
    mat4 model;
} pc;

layout(location = 0) in vec3 inPosition;
layout(location = 1) in vec3 inNormal;     // unused but keeps same vertex format
layout(location = 2) in float inQuality;   // unused

void main() {
    // Slight bias toward camera to prevent z-fighting with solid mesh
    vec4 pos = pc.mvp * vec4(inPosition, 1.0);
    pos.z -= 0.0001;
    gl_Position = pos;
}
