#include "mesh.h"
#include "vulkan_context.h"
#include <cstring>

VkVertexInputBindingDescription Vertex::bindingDescription() {
    VkVertexInputBindingDescription bd{};
    bd.binding = 0;
    bd.stride = sizeof(Vertex);
    bd.inputRate = VK_VERTEX_INPUT_RATE_VERTEX;
    return bd;
}

std::vector<VkVertexInputAttributeDescription> Vertex::attributeDescriptions() {
    return {
        {0, 0, VK_FORMAT_R32G32B32_SFLOAT, offsetof(Vertex, pos)},
        {1, 0, VK_FORMAT_R32G32B32_SFLOAT, offsetof(Vertex, normal)},
        {2, 0, VK_FORMAT_R32_SFLOAT,       offsetof(Vertex, quality)},
    };
}

static void createBuffer(VmaAllocator allocator, VkDeviceSize size,
                          VkBufferUsageFlags usage, VmaMemoryUsage memUsage,
                          VkBuffer& buffer, VmaAllocation& alloc) {
    VkBufferCreateInfo bci{};
    bci.sType = VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO;
    bci.size = size;
    bci.usage = usage;
    bci.sharingMode = VK_SHARING_MODE_EXCLUSIVE;

    VmaAllocationCreateInfo aci{};
    aci.usage = memUsage;

    vmaCreateBuffer(allocator, &bci, &aci, &buffer, &alloc, nullptr);
}

static void uploadToGPU(VulkanContext& ctx, const void* data, VkDeviceSize size,
                         VkBufferUsageFlags usage, VkBuffer& buffer, VmaAllocation& alloc) {
    // Create staging buffer
    VkBuffer staging;
    VmaAllocation stagingAlloc;
    createBuffer(ctx.allocator, size,
                 VK_BUFFER_USAGE_TRANSFER_SRC_BIT, VMA_MEMORY_USAGE_CPU_ONLY,
                 staging, stagingAlloc);

    void* mapped;
    vmaMapMemory(ctx.allocator, stagingAlloc, &mapped);
    memcpy(mapped, data, size);
    vmaUnmapMemory(ctx.allocator, stagingAlloc);

    // Create device-local buffer
    createBuffer(ctx.allocator, size,
                 usage | VK_BUFFER_USAGE_TRANSFER_DST_BIT, VMA_MEMORY_USAGE_GPU_ONLY,
                 buffer, alloc);

    // Copy
    VkCommandBuffer cmd = ctx.beginSingleTimeCommands();
    VkBufferCopy region{};
    region.size = size;
    vkCmdCopyBuffer(cmd, staging, buffer, 1, &region);
    ctx.endSingleTimeCommands(cmd);

    vmaDestroyBuffer(ctx.allocator, staging, stagingAlloc);
}

void GPUMesh::upload(VulkanContext& ctx,
                     const std::vector<Vertex>& vertices,
                     const std::vector<uint32_t>& triangleIndices,
                     const std::vector<uint32_t>& edgeIndices) {
    // Clean up old buffers if any
    cleanup(ctx.allocator);

    indexCount = (uint32_t)triangleIndices.size();
    edgeIndexCount = (uint32_t)edgeIndices.size();

    uploadToGPU(ctx, vertices.data(), vertices.size() * sizeof(Vertex),
                VK_BUFFER_USAGE_VERTEX_BUFFER_BIT, vertexBuffer, vertexAlloc);

    uploadToGPU(ctx, triangleIndices.data(), triangleIndices.size() * sizeof(uint32_t),
                VK_BUFFER_USAGE_INDEX_BUFFER_BIT, indexBuffer, indexAlloc);

    if (!edgeIndices.empty()) {
        uploadToGPU(ctx, edgeIndices.data(), edgeIndices.size() * sizeof(uint32_t),
                    VK_BUFFER_USAGE_INDEX_BUFFER_BIT, edgeIndexBuffer, edgeIndexAlloc);
    }
}

void GPUMesh::cleanup(VmaAllocator allocator) {
    if (vertexBuffer) { vmaDestroyBuffer(allocator, vertexBuffer, vertexAlloc); vertexBuffer = VK_NULL_HANDLE; }
    if (indexBuffer) { vmaDestroyBuffer(allocator, indexBuffer, indexAlloc); indexBuffer = VK_NULL_HANDLE; }
    if (edgeIndexBuffer) { vmaDestroyBuffer(allocator, edgeIndexBuffer, edgeIndexAlloc); edgeIndexBuffer = VK_NULL_HANDLE; }
    indexCount = 0;
    edgeIndexCount = 0;
}
