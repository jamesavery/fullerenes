#include "vulkan_context.h"
#include "renderer.h"
#include "mesh.h"
#include "camera.h"
#include "ui.h"
#include "fullerene_model.h"

#include <imgui.h>

#include <fullerenes/isomerdb.hh>
#include <fullerenes/buckygen-wrapper.hh>
#include <fullerenes/triangulation.hh>
#include <fullerenes/deltahedron.hh>
#include <fullerenes/buckinverse.hh>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <cstdio>
#include <cmath>
#include <vector>
#include <thread>
#include <mutex>
#include <atomic>
#include <chrono>

using namespace buckinverse;

// --- Snapshot: a frozen Deltahedron state for the step player ---

struct Snapshot {
    int step;
    std::string phase;
    std::vector<coord3d> points;
    std::vector<tri_t> triangles;
    std::vector<std::vector<int>> neighbours;
    int N;
    double gmax_L;
    double angle_relerr;
    int n_concave;
};

// --- Step Player: runs optimization in background, captures snapshots ---

struct StepPlayer {
    std::mutex mutex;
    std::vector<Snapshot> snapshots;
    std::atomic<bool> running{false};
    std::thread worker;

    void start(int fullereneN, int isomerIndex, PipelineMode mode) {
        if (running) return;
        {
            std::lock_guard<std::mutex> lock(mutex);
            snapshots.clear();
        }
        running = true;
        worker = std::thread([this, fullereneN, isomerIndex, mode]() {
            run(fullereneN, isomerIndex, mode);
        });
        worker.detach();
    }

    int numSnapshots() {
        std::lock_guard<std::mutex> lock(mutex);
        return (int)snapshots.size();
    }

    Snapshot getSnapshot(int idx) {
        std::lock_guard<std::mutex> lock(mutex);
        if (idx < 0 || idx >= (int)snapshots.size())
            return {};
        return snapshots[idx];
    }

private:
    void pushSnapshot(int step, const std::string& phase, const Deltahedron& D) {
        Snapshot snap;
        snap.step = step;
        snap.phase = phase;
        snap.points = D.points;
        snap.triangles = D.triangles;
        snap.neighbours = D.neighbours;
        snap.N = D.N;
        snap.gmax_L = D.final_gmax_L;
        snap.angle_relerr = D.final_angle_relerr;
        snap.n_concave = D.final_n_concave;
        std::lock_guard<std::mutex> lock(mutex);
        snapshots.push_back(std::move(snap));
    }

    // Find isomer and reduce to extension path (shared by all modes)
    bool findIsomer(int fullereneN, int isomerIndex, ExtensionPath& ep) {
        auto Q = BuckyGen::start(fullereneN, false);
        Graph G;
        int idx = 0;
        bool found = false;
        while (BuckyGen::next_fullerene(Q, G)) {
            if (idx == isomerIndex) { found = true; break; }
            idx++;
        }
        BuckyGen::stop(Q);
        if (!found) {
            fprintf(stderr, "Isomer C%d #%d not found\n", fullereneN, isomerIndex);
            return false;
        }
        ReducibleDual rd(G);
        ep = rd.reduceToExtensionPath(20);
        fprintf(stderr, "C%d idx=%d seed=%d steps=%zu\n",
                fullereneN, isomerIndex, (int)ep.seed, ep.steps.size());
        return true;
    }

    void runRaw(const ExtensionPath& ep) {
        // Build unoptimized geometry using fromExtensionPath on partial paths
        for (int k = 0; k <= (int)ep.steps.size(); k++) {
            ExtensionPath partial = ep;
            partial.steps.resize(k);
            Deltahedron D = Deltahedron::fromExtensionPath(partial);
            std::string label = (k == 0) ? "seed" : "step " + std::to_string(k);
            pushSnapshot(k, label, D);
        }
    }

    void runOptimized(const ExtensionPath& ep) {
        // fromExtensionPathOptimized with step callback (one snapshot per phase)
        auto callback = [this](int step, const char* phase, const Deltahedron& D) {
            pushSnapshot(step, phase, D);
        };
        Deltahedron D = Deltahedron::fromExtensionPathOptimized(
            ep, stderr, callback,
            OptMethod::LBFGS, 1e-4, 1e-11, 0, 0, 0,
            OptMethod::STEIHAUG
        );
        pushSnapshot(-1, "final", D);
    }

    void runLive(const ExtensionPath& ep) {
        // Phase-boundary snapshots via StepCallback
        auto stepCb = [this](int step, const char* phase, const Deltahedron& D) {
            pushSnapshot(step, phase, D);
        };

        // Per-iteration snapshots via IterCallback (propagated through the library)
        auto iterCb = [this](int iter, double gmax_L, double angle_relerr,
                             int n_concave, const Deltahedron& D) {
            Snapshot snap;
            snap.step = -1;
            snap.phase = "iter " + std::to_string(iter);
            snap.points = D.points;
            snap.triangles = D.triangles;
            snap.neighbours = D.neighbours;
            snap.N = D.N;
            snap.gmax_L = gmax_L;
            snap.angle_relerr = angle_relerr;
            snap.n_concave = n_concave;
            std::lock_guard<std::mutex> lock(mutex);
            snapshots.push_back(std::move(snap));
        };

        Deltahedron D = Deltahedron::fromExtensionPathOptimized(
            ep, stderr, stepCb,
            OptMethod::LBFGS, 1e-4, 1e-11, 0, 0, 0,
            OptMethod::STEIHAUG, 1e-10, false,
            iterCb, 5  // every 5 iterations
        );
        pushSnapshot(-1, "final", D);
    }

    void run(int fullereneN, int isomerIndex, PipelineMode mode) {
        ExtensionPath ep;
        if (!findIsomer(fullereneN, isomerIndex, ep)) {
            running = false;
            return;
        }

        switch (mode) {
            case PipelineMode::Raw:       runRaw(ep); break;
            case PipelineMode::Optimized: runOptimized(ep); break;
            case PipelineMode::Live:      runLive(ep); break;
        }

        running = false;
        fprintf(stderr, "Done. %d snapshots captured.\n", numSnapshots());
    }
};

// --- Build a temporary Deltahedron from a snapshot (for mesh conversion) ---

static Deltahedron snapshotToDeltahedron(const Snapshot& snap) {
    Triangulation T;
    T.N = snap.N;
    T.neighbours = snap.neighbours;
    T.is_oriented = true;
    T.triangles = snap.triangles;
    Deltahedron D(T, snap.points);
    D.final_gmax_L = snap.gmax_L;
    D.final_angle_relerr = snap.angle_relerr;
    D.final_n_concave = snap.n_concave;
    return D;
}

// --- Compute quality stats from a Deltahedron ---

static void computeQualityStats(const Deltahedron& D, UIState& state) {
    state.vertexCount = D.N;
    state.triangleCount = (int)D.triangles.size();
    state.gmaxL = (float)D.final_gmax_L;
    state.concaveCount = D.final_n_concave;

    // Edge CV
    double sum = 0, sum2 = 0; int ne = 0;
    for (int u = 0; u < D.N; u++)
        for (int v : D.neighbours[u])
            if (v > u) {
                double l = (D.points[u] - D.points[v]).norm();
                sum += l; sum2 += l*l; ne++;
            }
    if (ne > 0) {
        double L = sum / ne;
        state.edgeCv = (float)(sqrt(std::max(0.0, sum2/ne - L*L)) / L);
    }

    // Angle range
    double amin = 180, amax = 0;
    for (const auto& tri : D.triangles) {
        for (int c = 0; c < 3; c++) {
            coord3d va = D.points[tri[(c+1)%3]] - D.points[tri[c]];
            coord3d vb = D.points[tri[(c+2)%3]] - D.points[tri[c]];
            double ang = coord3d::angle(va, vb) * 180.0 / M_PI;
            amin = std::min(amin, ang);
            amax = std::max(amax, ang);
        }
    }
    state.angleMin = (float)amin;
    state.angleMax = (float)amax;
}

// --- Main ---

int main(int argc, char** argv) {
    VulkanContext ctx;
    ctx.init(1280, 800, "Fullerene Viewer");

    Renderer renderer;
    renderer.init(ctx);

    UI ui;
    ui.init(ctx, renderer);

    Camera camera;
    // Input state for polling (no GLFW callbacks — ImGui owns those)
    bool prevLMB = false, prevRMB = false;
    double prevMouseX = 0, prevMouseY = 0;
    double prevScrollY = 0;

    GPUMesh mesh;
    UIState uiState;
    StepPlayer player;
    bool meshLoaded = false;
    bool needCameraFit = true;
    int lastSnapshotIndex = -1;
    ColorMode lastColorMode = ColorMode::Solid;

    // Auto-load C60 IPR isomer (buckminsterfullerene) at startup
    uiState.fullereneN = 60;
    uiState.isomerIndex = 0;
    uiState.loadRequested = true;

    auto lastFrameTime = std::chrono::steady_clock::now();
    float stepAccumulator = 0;

    while (!glfwWindowShouldClose(ctx.window)) {
        glfwPollEvents();

        // --- Poll mouse input, forwarding to camera only when ImGui doesn't want it ---
        {
            ImGuiIO& io = ImGui::GetIO();
            double mx, my;
            glfwGetCursorPos(ctx.window, &mx, &my);

            if (!io.WantCaptureMouse) {
                bool lmb = glfwGetMouseButton(ctx.window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
                bool rmb = glfwGetMouseButton(ctx.window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS;

                // Detect button press/release edges
                if (lmb && !prevLMB) camera.handleMouseButton(0, 1, mx, my);
                if (!lmb && prevLMB) camera.handleMouseButton(0, 0, mx, my);
                if (rmb && !prevRMB) camera.handleMouseButton(1, 1, mx, my);
                if (!rmb && prevRMB) camera.handleMouseButton(1, 0, mx, my);

                // Mouse move
                if (mx != prevMouseX || my != prevMouseY) {
                    camera.handleMouseMove(mx, my);
                }

                // Scroll (read from ImGui's accumulated scroll)
                if (io.MouseWheel != 0) {
                    camera.handleScroll((double)io.MouseWheel);
                }

                prevLMB = lmb;
                prevRMB = rmb;
            } else {
                // ImGui has focus — release any camera drag
                if (prevLMB) { camera.handleMouseButton(0, 0, mx, my); prevLMB = false; }
                if (prevRMB) { camera.handleMouseButton(1, 0, mx, my); prevRMB = false; }
            }
            prevMouseX = mx;
            prevMouseY = my;
        }

        auto now = std::chrono::steady_clock::now();
        float dt = std::chrono::duration<float>(now - lastFrameTime).count();
        lastFrameTime = now;

        // Handle load request
        if (uiState.loadRequested) {
            uiState.loadRequested = false;
            player.start(uiState.fullereneN, uiState.isomerIndex, uiState.pipelineMode);
            uiState.currentStep = 0;
            uiState.totalSteps = 0;
            meshLoaded = false;
            needCameraFit = true;
            lastSnapshotIndex = -1;
        }

        // Update step count from player
        int numSnaps = player.numSnapshots();
        if (numSnaps > 0) {
            uiState.totalSteps = numSnaps - 1;
        }

        // Auto-play
        if (uiState.playing && uiState.totalSteps > 0) {
            stepAccumulator += dt * uiState.playSpeed;
            while (stepAccumulator >= 1.0f) {
                stepAccumulator -= 1.0f;
                uiState.currentStep++;
                if (uiState.currentStep > uiState.totalSteps) {
                    uiState.currentStep = uiState.totalSteps;
                    uiState.playing = false;
                }
            }
        }

        // Upload new mesh if snapshot changed or color mode changed
        int targetSnap = std::min(uiState.currentStep, numSnaps - 1);
        if (numSnaps > 0 && (targetSnap != lastSnapshotIndex || uiState.colorMode != lastColorMode)) {
            lastSnapshotIndex = targetSnap;
            lastColorMode = uiState.colorMode;

            Snapshot snap = player.getSnapshot(targetSnap);
            if (snap.N > 0) {
                Deltahedron D = snapshotToDeltahedron(snap);
                MeshData data = deltahedronToMesh(D, uiState.colorMode);
                mesh.upload(ctx, data.vertices, data.triangleIndices, data.edgeIndices);
                meshLoaded = true;

                // Always keep orbit center at mesh barycentre
                camera.target = data.center;

                // Auto-fit camera distance on first load
                if (needCameraFit) {
                    camera.fitToBounds(data.center, data.radius);
                    needCameraFit = false;
                }

                computeQualityStats(D, uiState);

                // Symmetry check (skip for mid-iteration snapshots — too expensive)
                if (snap.phase.substr(0, 4) != "iter") {
                    SymmetryCheck sc = checkSymmetry(D);
                    uiState.pointGroup = sc.pointGroup;
                    uiState.symmetryOrder = sc.groupOrder;
                    uiState.geometricSymCount = sc.geometricCount;
                    uiState.symmetryMaxRMSD = sc.maxRMSD;
                    uiState.symmetryMeanRMSD = sc.meanRMSD;
                }

                // Build snapshot label
                char label[128];
                if (snap.step >= 0) {
                    snprintf(label, sizeof(label), "Step %d: %s (N=%d)",
                             snap.step, snap.phase.c_str(), snap.N);
                } else {
                    snprintf(label, sizeof(label), "%s (N=%d)", snap.phase.c_str(), snap.N);
                }
                uiState.snapshotLabel = label;
            }
        }

        // --- Render frame ---

        // Handle minimized window
        int fw, fh;
        glfwGetFramebufferSize(ctx.window, &fw, &fh);
        if (fw == 0 || fh == 0) continue;

        // Wait for previous frame
        vkWaitForFences(ctx.device, 1, &ctx.inFlightFence[ctx.currentFrame], VK_TRUE, UINT64_MAX);

        uint32_t imageIndex;
        VkResult result = vkAcquireNextImageKHR(ctx.device, ctx.swapchain, UINT64_MAX,
                                                 ctx.imageAvailableSem[ctx.currentFrame],
                                                 VK_NULL_HANDLE, &imageIndex);
        if (result == VK_ERROR_OUT_OF_DATE_KHR) {
            ctx.recreateSwapchain();
            renderer.recreateFramebuffers(ctx);
            continue;
        }

        vkResetFences(ctx.device, 1, &ctx.inFlightFence[ctx.currentFrame]);

        VkCommandBuffer cmd = ctx.commandBuffers[ctx.currentFrame];
        vkResetCommandBuffer(cmd, 0);

        // Update UBO
        float aspect = (float)ctx.swapchainExtent.width / (float)ctx.swapchainExtent.height;
        SceneUBO ubo{};
        ubo.lightDir = glm::normalize(glm::vec3(0.5f, 1.0f, 0.8f));
        ubo.eyePos = camera.eyePosition();
        ubo.baseColor = glm::vec4(0.4f, 0.65f, 0.85f, 1.0f);  // steel blue
        ubo.colorMode = (int)uiState.colorMode;
        renderer.updateUBO(ctx.currentFrame, ubo);

        // Begin render pass
        renderer.beginFrame(ctx, cmd, imageIndex);

        if (meshLoaded) {
            // Compute MVP
            glm::mat4 view = camera.viewMatrix();
            glm::mat4 proj = camera.projectionMatrix(aspect);
            glm::mat4 model = glm::mat4(1.0f);

            PushConstants pc{};
            pc.mvp = proj * view * model;
            pc.model = model;

            vkCmdPushConstants(cmd, renderer.pipelineLayout, VK_SHADER_STAGE_VERTEX_BIT,
                               0, sizeof(PushConstants), &pc);

            if (uiState.showSolid) {
                renderer.drawMesh(cmd, mesh, false, ctx.currentFrame);
            }
            if (uiState.showWireframe) {
                renderer.drawMesh(cmd, mesh, true, ctx.currentFrame);
            }
        }

        // ImGui
        ui.newFrame();
        ui.draw(uiState);
        ui.render(cmd);

        renderer.endFrame(cmd);

        // Submit
        VkSemaphore waitSems[] = {ctx.imageAvailableSem[ctx.currentFrame]};
        VkPipelineStageFlags waitStages[] = {VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT};
        VkSemaphore signalSems[] = {ctx.renderFinishedSem[ctx.currentFrame]};

        VkSubmitInfo submitInfo{};
        submitInfo.sType = VK_STRUCTURE_TYPE_SUBMIT_INFO;
        submitInfo.waitSemaphoreCount = 1;
        submitInfo.pWaitSemaphores = waitSems;
        submitInfo.pWaitDstStageMask = waitStages;
        submitInfo.commandBufferCount = 1;
        submitInfo.pCommandBuffers = &cmd;
        submitInfo.signalSemaphoreCount = 1;
        submitInfo.pSignalSemaphores = signalSems;

        if (vkQueueSubmit(ctx.graphicsQueue, 1, &submitInfo, ctx.inFlightFence[ctx.currentFrame]) != VK_SUCCESS)
            throw std::runtime_error("Failed to submit draw command buffer");

        VkPresentInfoKHR presentInfo{};
        presentInfo.sType = VK_STRUCTURE_TYPE_PRESENT_INFO_KHR;
        presentInfo.waitSemaphoreCount = 1;
        presentInfo.pWaitSemaphores = signalSems;
        presentInfo.swapchainCount = 1;
        presentInfo.pSwapchains = &ctx.swapchain;
        presentInfo.pImageIndices = &imageIndex;

        result = vkQueuePresentKHR(ctx.presentQueue, &presentInfo);
        if (result == VK_ERROR_OUT_OF_DATE_KHR || result == VK_SUBOPTIMAL_KHR || ctx.framebufferResized) {
            ctx.recreateSwapchain();
            renderer.recreateFramebuffers(ctx);
        }

        ctx.currentFrame = (ctx.currentFrame + 1) % VulkanContext::MAX_FRAMES_IN_FLIGHT;
    }

    vkDeviceWaitIdle(ctx.device);
    mesh.cleanup(ctx.allocator);
    ui.cleanup(ctx);
    renderer.cleanup(ctx);
    ctx.cleanup();
    return 0;
}
