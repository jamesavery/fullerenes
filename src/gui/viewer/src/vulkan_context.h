#pragma once

#include <vulkan/vulkan.h>
#include <vk_mem_alloc.h>
#include <GLFW/glfw3.h>
#include <vector>
#include <stdexcept>
#include <functional>

struct VulkanContext {
    GLFWwindow* window = nullptr;
    VkInstance instance = VK_NULL_HANDLE;
    VkSurfaceKHR surface = VK_NULL_HANDLE;
    VkPhysicalDevice physicalDevice = VK_NULL_HANDLE;
    VkDevice device = VK_NULL_HANDLE;
    VkQueue graphicsQueue = VK_NULL_HANDLE;
    VkQueue presentQueue = VK_NULL_HANDLE;
    uint32_t graphicsFamily = 0;
    uint32_t presentFamily = 0;
    VmaAllocator allocator = VK_NULL_HANDLE;

    // Swapchain
    VkSwapchainKHR swapchain = VK_NULL_HANDLE;
    VkFormat swapchainFormat = VK_FORMAT_UNDEFINED;
    VkExtent2D swapchainExtent = {};
    std::vector<VkImage> swapchainImages;
    std::vector<VkImageView> swapchainImageViews;

    // Depth buffer
    VkImage depthImage = VK_NULL_HANDLE;
    VmaAllocation depthAlloc = VK_NULL_HANDLE;
    VkImageView depthImageView = VK_NULL_HANDLE;
    VkFormat depthFormat = VK_FORMAT_D32_SFLOAT;

    // Command pool
    VkCommandPool commandPool = VK_NULL_HANDLE;

    // Sync
    static constexpr int MAX_FRAMES_IN_FLIGHT = 2;
    VkSemaphore imageAvailableSem[MAX_FRAMES_IN_FLIGHT] = {};
    VkSemaphore renderFinishedSem[MAX_FRAMES_IN_FLIGHT] = {};
    VkFence inFlightFence[MAX_FRAMES_IN_FLIGHT] = {};
    VkCommandBuffer commandBuffers[MAX_FRAMES_IN_FLIGHT] = {};
    int currentFrame = 0;

    bool framebufferResized = false;

    void init(int width, int height, const char* title);
    void cleanup();
    void recreateSwapchain();

    // One-shot command buffer helper
    VkCommandBuffer beginSingleTimeCommands();
    void endSingleTimeCommands(VkCommandBuffer cmd);

private:
    void createInstance();
    void createSurface();
    void pickPhysicalDevice();
    void createDevice();
    void createAllocator();
    void createSwapchain();
    void createDepthResources();
    void createCommandPool();
    void createSyncObjects();
    void cleanupSwapchain();
};
