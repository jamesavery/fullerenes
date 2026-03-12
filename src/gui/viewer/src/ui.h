#pragma once

#include <vulkan/vulkan.h>
#include <string>
#include "fullerene_model.h"

struct VulkanContext;
struct Renderer;

// Pipeline mode: what to show as snapshots
enum class PipelineMode {
    Raw = 0,       // Unoptimized: just strip placement (fromExtensionPath)
    Optimized = 1, // Per-step optimized (fromExtensionPathOptimized), one snapshot per phase
    Live = 2,      // Per-step optimized + per-iteration snapshots during convergence
};

struct UIState {
    int fullereneN = 60;
    int isomerIndex = 0;
    bool loadRequested = false;
    bool showWireframe = true;
    bool showSolid = true;
    ColorMode colorMode = ColorMode::Solid;
    PipelineMode pipelineMode = PipelineMode::Live;

    // Step playback
    int currentStep = 0;
    int totalSteps = 0;
    bool playing = false;
    float playSpeed = 2.0f;  // steps per second

    // Quality stats display
    float edgeCv = 0;
    float angleMin = 0, angleMax = 0;
    int concaveCount = 0;
    float gmaxL = 0;
    int vertexCount = 0;
    int triangleCount = 0;

    // Symmetry validation
    std::string pointGroup;
    int symmetryOrder = 0;
    int geometricSymCount = 0;
    double symmetryMaxRMSD = 0;
    double symmetryMeanRMSD = 0;

    // Current snapshot label
    std::string snapshotLabel;
};

struct UI {
    VkDescriptorPool imguiPool = VK_NULL_HANDLE;

    void init(VulkanContext& ctx, Renderer& renderer);
    void cleanup(VulkanContext& ctx);
    void newFrame();
    void render(VkCommandBuffer cmd);
    void draw(UIState& state);
};
