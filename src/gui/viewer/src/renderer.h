#pragma once

#include <vulkan/vulkan.h>
#include <vk_mem_alloc.h>
#include <glm/glm.hpp>
#include <vector>
#include <string>

struct VulkanContext;
struct GPUMesh;

struct PushConstants {
    glm::mat4 mvp;
    glm::mat4 model;
};

struct SceneUBO {
    glm::vec3 lightDir;
    float _pad0;
    glm::vec3 eyePos;
    float _pad1;
    glm::vec4 baseColor;
    int colorMode;
    float _pad2[3];
};

struct Renderer {
    VkRenderPass renderPass = VK_NULL_HANDLE;
    VkPipelineLayout pipelineLayout = VK_NULL_HANDLE;
    VkPipeline solidPipeline = VK_NULL_HANDLE;
    VkPipeline wireframePipeline = VK_NULL_HANDLE;

    VkDescriptorSetLayout descriptorSetLayout = VK_NULL_HANDLE;
    VkDescriptorPool descriptorPool = VK_NULL_HANDLE;
    std::vector<VkDescriptorSet> descriptorSets;

    // Per-frame UBO
    std::vector<VkBuffer> uboBuffers;
    std::vector<VmaAllocation> uboAllocs;
    std::vector<void*> uboMapped;

    std::vector<VkFramebuffer> framebuffers;

    void init(VulkanContext& ctx);
    void cleanup(VulkanContext& ctx);
    void recreateFramebuffers(VulkanContext& ctx);

    void beginFrame(VulkanContext& ctx, VkCommandBuffer cmd, uint32_t imageIndex);
    void endFrame(VkCommandBuffer cmd);

    void drawMesh(VkCommandBuffer cmd, const GPUMesh& mesh, bool wireframe, int frameIndex);
    void updateUBO(int frameIndex, const SceneUBO& ubo);

private:
    VkShaderModule loadShaderModule(VkDevice device, const std::string& path);
    void createRenderPass(VulkanContext& ctx);
    void createDescriptorSetLayout(VkDevice device);
    void createPipelines(VulkanContext& ctx);
    void createUBOs(VulkanContext& ctx);
    void createDescriptorSets(VulkanContext& ctx);
};
