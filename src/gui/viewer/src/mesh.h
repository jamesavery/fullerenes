#pragma once

#include <vulkan/vulkan.h>
#include <vk_mem_alloc.h>
#include <glm/glm.hpp>
#include <vector>

struct VulkanContext;

struct Vertex {
    glm::vec3 pos;
    glm::vec3 normal;
    float quality;  // per-vertex metric for color-coding

    static VkVertexInputBindingDescription bindingDescription();
    static std::vector<VkVertexInputAttributeDescription> attributeDescriptions();
};

struct GPUMesh {
    VkBuffer vertexBuffer = VK_NULL_HANDLE;
    VmaAllocation vertexAlloc = VK_NULL_HANDLE;
    VkBuffer indexBuffer = VK_NULL_HANDLE;
    VmaAllocation indexAlloc = VK_NULL_HANDLE;
    uint32_t indexCount = 0;

    // Wireframe edge indices (line list)
    VkBuffer edgeIndexBuffer = VK_NULL_HANDLE;
    VmaAllocation edgeIndexAlloc = VK_NULL_HANDLE;
    uint32_t edgeIndexCount = 0;

    void upload(VulkanContext& ctx,
                const std::vector<Vertex>& vertices,
                const std::vector<uint32_t>& triangleIndices,
                const std::vector<uint32_t>& edgeIndices);
    void cleanup(VmaAllocator allocator);
};
