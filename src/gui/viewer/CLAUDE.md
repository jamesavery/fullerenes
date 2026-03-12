# Fullerene Viewer - CLAUDE.md

## Overview

Vulkan-accelerated real-time viewer for fullerene deltahedra, integrated with the extension path optimizer. Displays 3D triangle meshes with Blinn-Phong shading, wireframe overlay, and quality heat maps.

## Build

```bash
cd viewer
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build .
```

Requires: Vulkan SDK (MoltenVK on macOS), GLFW, GLM, all via Homebrew. ImGui and VMA are fetched automatically via CMake FetchContent.

The fullerene library is found relative to this directory (../../build/src/c++/libfullerenes.dylib). Override with `-DFULLERENE_ROOT=... -DFULLERENE_BUILD_DIR=...`.

## Run

```bash
FULLERENE_DATABASE_PATH=/path/to/fullerenes/database ./viewer
```

## Architecture

- `vulkan_context.h/.cpp` -- VkInstance, VkDevice, swapchain, VMA allocator, sync primitives
- `renderer.h/.cpp` -- Render pass, solid + wireframe pipelines, descriptor sets, UBOs
- `mesh.h/.cpp` -- GPU vertex/index buffer upload via VMA staging
- `camera.h/.cpp` -- Orbit camera (azimuth/elevation/distance), pan, zoom, auto-fit
- `ui.h/.cpp` -- Dear ImGui overlay: isomer selector, display toggles, step playback, quality stats
- `fullerene_model.h/.cpp` -- Bridge: Deltahedron -> GPU vertex/index data with normals and quality metrics
- `main.cpp` -- Entry point, main loop, StepPlayer (background optimizer thread with snapshot capture)

## Shaders

- `mesh.vert/frag` -- MVP transform + Blinn-Phong with heat map coloring
- `wireframe.vert/frag` -- Edge line overlay with depth bias

## Key Design Decisions

- StepPlayer runs optimization in a background thread; the main thread renders snapshots captured via the Deltahedron::StepCallback mechanism
- Vertex format: position (vec3) + normal (vec3) + quality (float) -- normals are area-weighted face normals averaged to vertices
- Edge wireframe uses LINE_LIST topology with explicit edge index buffer built from the adjacency list
- Double-to-float conversion happens in fullerene_model.cpp (fullerene geometry is O(1) scale, float is sufficient)
- MoltenVK portability: VK_INSTANCE_CREATE_ENUMERATE_PORTABILITY_BIT_KHR and VK_KHR_portability_subset device extension are enabled

## Color Modes

0. Solid color (steel blue)
1. Angle error -- max |theta - 60| / 60 per vertex, normalized to [0,1], heat mapped
2. Convexity -- signed height h, concave (h < 0) maps to red
3. Degree -- pentagon vertices (degree 5) highlighted
