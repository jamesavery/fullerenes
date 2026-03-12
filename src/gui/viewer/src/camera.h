#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

struct Camera {
    float distance = 5.0f;
    float azimuth = 0.0f;     // radians
    float elevation = 0.3f;   // radians
    glm::vec3 target = {0, 0, 0};
    float fov = 45.0f;
    float nearClip = 0.01f;
    float farClip = 100.0f;

    // Mouse state
    bool rotating = false;
    bool panning = false;
    double lastX = 0, lastY = 0;

    glm::vec3 eyePosition() const;
    glm::mat4 viewMatrix() const;
    glm::mat4 projectionMatrix(float aspect) const;

    void handleMouseButton(int button, int action, double x, double y);
    void handleMouseMove(double x, double y);
    void handleScroll(double yoffset);

    // Auto-fit camera to mesh bounding sphere
    void fitToBounds(glm::vec3 center, float radius);
};
