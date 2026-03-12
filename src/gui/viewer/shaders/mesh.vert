#version 450

layout(push_constant) uniform PushConstants {
    mat4 mvp;
    mat4 model;  // for transforming normals
} pc;

layout(location = 0) in vec3 inPosition;
layout(location = 1) in vec3 inNormal;
layout(location = 2) in float inQuality;  // per-vertex quality metric

layout(location = 0) out vec3 fragNormal;
layout(location = 1) out vec3 fragWorldPos;
layout(location = 2) out float fragQuality;

void main() {
    gl_Position = pc.mvp * vec4(inPosition, 1.0);
    fragWorldPos = (pc.model * vec4(inPosition, 1.0)).xyz;
    fragNormal = mat3(pc.model) * inNormal;
    fragQuality = inQuality;
}
