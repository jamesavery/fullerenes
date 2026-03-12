#include "ui.h"
#include "vulkan_context.h"
#include "renderer.h"
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_vulkan.h>

void UI::init(VulkanContext& ctx, Renderer& renderer) {
    // Create dedicated descriptor pool for ImGui
    VkDescriptorPoolSize poolSizes[] = {
        {VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1},
    };
    VkDescriptorPoolCreateInfo dpci{};
    dpci.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO;
    dpci.flags = VK_DESCRIPTOR_POOL_CREATE_FREE_DESCRIPTOR_SET_BIT;
    dpci.maxSets = 1;
    dpci.poolSizeCount = 1;
    dpci.pPoolSizes = poolSizes;
    vkCreateDescriptorPool(ctx.device, &dpci, nullptr, &imguiPool);

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGui::StyleColorsDark();

    ImGui_ImplGlfw_InitForVulkan(ctx.window, true);

    ImGui_ImplVulkan_InitInfo initInfo{};
    initInfo.Instance = ctx.instance;
    initInfo.PhysicalDevice = ctx.physicalDevice;
    initInfo.Device = ctx.device;
    initInfo.QueueFamily = ctx.graphicsFamily;
    initInfo.Queue = ctx.graphicsQueue;
    initInfo.DescriptorPool = imguiPool;
    initInfo.MinImageCount = 2;
    initInfo.ImageCount = (uint32_t)ctx.swapchainImages.size();
    initInfo.RenderPass = renderer.renderPass;
    initInfo.MSAASamples = VK_SAMPLE_COUNT_1_BIT;

    ImGui_ImplVulkan_Init(&initInfo);

    // Upload fonts
    ImGui_ImplVulkan_CreateFontsTexture();
}

void UI::cleanup(VulkanContext& ctx) {
    ImGui_ImplVulkan_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    vkDestroyDescriptorPool(ctx.device, imguiPool, nullptr);
}

void UI::newFrame() {
    ImGui_ImplVulkan_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
}

void UI::render(VkCommandBuffer cmd) {
    ImGui::Render();
    ImGui_ImplVulkan_RenderDrawData(ImGui::GetDrawData(), cmd);
}

void UI::draw(UIState& state) {
    ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(300, 400), ImGuiCond_FirstUseEver);

    ImGui::Begin("Fullerene Viewer");

    // Isomer loading
    ImGui::SeparatorText("Isomer");
    ImGui::InputInt("N (vertices)", &state.fullereneN, 2, 10);
    if (state.fullereneN < 20) state.fullereneN = 20;
    if (state.fullereneN % 2 != 0) state.fullereneN++;
    ImGui::InputInt("Index", &state.isomerIndex);
    if (state.isomerIndex < 0) state.isomerIndex = 0;

    // Pipeline mode
    const char* pipelineModes[] = {"Raw (no optimization)", "Optimized (per-step)", "Live (per-iteration)"};
    int pm = (int)state.pipelineMode;
    ImGui::Combo("Mode", &pm, pipelineModes, 3);
    state.pipelineMode = (PipelineMode)pm;

    if (ImGui::Button("Load Isomer")) {
        state.loadRequested = true;
    }

    // Rendering options
    ImGui::SeparatorText("Display");
    ImGui::Checkbox("Solid", &state.showSolid);
    ImGui::SameLine();
    ImGui::Checkbox("Wireframe", &state.showWireframe);

    const char* colorModes[] = {"Solid Color", "Angle Error", "Convexity", "Degree (5 vs 6)"};
    int cm = (int)state.colorMode;
    if (ImGui::Combo("Color", &cm, colorModes, 4)) {
        state.colorMode = (ColorMode)cm;
    }

    // Step playback
    if (state.totalSteps > 0) {
        ImGui::SeparatorText("Steps");
        ImGui::SliderInt("Step", &state.currentStep, 0, state.totalSteps);
        if (!state.snapshotLabel.empty()) {
            ImGui::TextWrapped("%s", state.snapshotLabel.c_str());
        }
        if (ImGui::Button(state.playing ? "Pause" : "Play")) {
            state.playing = !state.playing;
        }
        ImGui::SameLine();
        if (ImGui::Button("|<")) state.currentStep = 0;
        ImGui::SameLine();
        if (ImGui::Button(">|")) state.currentStep = state.totalSteps;
        ImGui::SliderFloat("Speed", &state.playSpeed, 0.5f, 30.0f);
    }

    // Quality stats
    ImGui::SeparatorText("Quality");
    ImGui::Text("Vertices: %d", state.vertexCount);
    ImGui::Text("Triangles: %d", state.triangleCount);
    ImGui::Text("Edge CV: %.3e", state.edgeCv);
    ImGui::Text("Angle: [%.2f, %.2f] deg", state.angleMin, state.angleMax);
    ImGui::Text("Concave: %d", state.concaveCount);
    ImGui::Text("gmax*L: %.2e", state.gmaxL);

    // Symmetry validation
    if (state.symmetryOrder > 0) {
        ImGui::SeparatorText("Symmetry");
        ImGui::Text("Point group: %s (|G|=%d)", state.pointGroup.c_str(), state.symmetryOrder);
        bool allMatch = (state.geometricSymCount == state.symmetryOrder);
        if (allMatch) {
            ImGui::TextColored(ImVec4(0.3f, 1.0f, 0.3f, 1.0f), "Geometric: %d/%d",
                               state.geometricSymCount, state.symmetryOrder);
        } else {
            ImGui::TextColored(ImVec4(1.0f, 0.3f, 0.3f, 1.0f), "Geometric: %d/%d BROKEN",
                               state.geometricSymCount, state.symmetryOrder);
        }
        ImGui::Text("Max RMSD: %.2e", state.symmetryMaxRMSD);
        ImGui::Text("Mean RMSD: %.2e", state.symmetryMeanRMSD);
    }

    ImGui::End();
}
