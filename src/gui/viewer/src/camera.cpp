#include "camera.h"
#include <algorithm>
#include <cmath>

glm::vec3 Camera::eyePosition() const {
    float x = distance * cos(elevation) * sin(azimuth);
    float y = distance * sin(elevation);
    float z = distance * cos(elevation) * cos(azimuth);
    return target + glm::vec3(x, y, z);
}

glm::mat4 Camera::viewMatrix() const {
    return glm::lookAt(eyePosition(), target, glm::vec3(0, 1, 0));
}

glm::mat4 Camera::projectionMatrix(float aspect) const {
    return glm::perspective(glm::radians(fov), aspect, nearClip, farClip);
}

void Camera::handleMouseButton(int button, int action, double x, double y) {
    if (button == 0) { // left
        rotating = (action == 1);
        lastX = x; lastY = y;
    } else if (button == 1) { // right
        panning = (action == 1);
        lastX = x; lastY = y;
    }
}

void Camera::handleMouseMove(double x, double y) {
    double dx = x - lastX;
    double dy = y - lastY;
    lastX = x; lastY = y;

    if (rotating) {
        azimuth -= (float)dx * 0.005f;
        elevation += (float)dy * 0.005f;
        elevation = std::clamp(elevation, -1.5f, 1.5f);
    }
    if (panning) {
        // Pan in camera-local XY plane
        glm::vec3 eye = eyePosition();
        glm::vec3 forward = glm::normalize(target - eye);
        glm::vec3 right = glm::normalize(glm::cross(forward, glm::vec3(0, 1, 0)));
        glm::vec3 up = glm::cross(right, forward);

        float scale = distance * 0.002f;
        target -= right * (float)dx * scale;
        target += up * (float)dy * scale;
    }
}

void Camera::handleScroll(double yoffset) {
    distance *= (1.0f - (float)yoffset * 0.1f);
    distance = std::max(0.1f, distance);
}

void Camera::fitToBounds(glm::vec3 center, float radius) {
    target = center;
    distance = radius * 2.5f / tan(glm::radians(fov * 0.5f));
    nearClip = std::max(0.001f, distance - radius * 3.0f);
    farClip = distance + radius * 3.0f;
}
