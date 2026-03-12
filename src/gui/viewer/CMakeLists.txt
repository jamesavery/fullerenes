# Vulkan-accelerated fullerene/deltahedron viewer
# Optional: only built when Vulkan, GLFW, GLM, and LAPACK are available.

find_package(Vulkan QUIET)
find_package(glfw3 QUIET)
find_package(glm QUIET)
find_package(LAPACK QUIET)

if(NOT Vulkan_FOUND OR NOT glfw3_FOUND OR NOT glm_FOUND OR NOT LAPACK_FOUND)
    message(STATUS "Viewer: skipping (missing Vulkan/GLFW/GLM/LAPACK)")
    return()
endif()

message(STATUS "Viewer: building fullerene-viewer")

# Dear ImGui via FetchContent
include(FetchContent)

FetchContent_Declare(
    imgui
    GIT_REPOSITORY https://github.com/ocornut/imgui.git
    GIT_TAG        v1.91.8
)
FetchContent_MakeAvailable(imgui)

add_library(imgui STATIC
    ${imgui_SOURCE_DIR}/imgui.cpp
    ${imgui_SOURCE_DIR}/imgui_demo.cpp
    ${imgui_SOURCE_DIR}/imgui_draw.cpp
    ${imgui_SOURCE_DIR}/imgui_tables.cpp
    ${imgui_SOURCE_DIR}/imgui_widgets.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_vulkan.cpp
)
target_include_directories(imgui PUBLIC
    ${imgui_SOURCE_DIR}
    ${imgui_SOURCE_DIR}/backends
)
target_link_libraries(imgui PUBLIC Vulkan::Vulkan glfw)

# VMA (header-only)
FetchContent_Declare(
    vma
    GIT_REPOSITORY https://github.com/GPUOpen-LibrariesAndSDKs/VulkanMemoryAllocator.git
    GIT_TAG        v3.2.1
)
FetchContent_MakeAvailable(vma)

# --- Shader compilation ---

find_program(GLSLC glslc HINTS /opt/homebrew/bin)
if(NOT GLSLC)
    message(STATUS "Viewer: skipping (glslc not found)")
    return()
endif()

set(VIEWER_SHADER_DIR ${CMAKE_CURRENT_SOURCE_DIR}/shaders)
set(VIEWER_SPIRV_DIR ${CMAKE_CURRENT_BINARY_DIR}/shaders)
file(MAKE_DIRECTORY ${VIEWER_SPIRV_DIR})

set(VIEWER_SHADERS
    ${VIEWER_SHADER_DIR}/mesh.vert
    ${VIEWER_SHADER_DIR}/mesh.frag
    ${VIEWER_SHADER_DIR}/wireframe.vert
    ${VIEWER_SHADER_DIR}/wireframe.frag
)

foreach(SHADER ${VIEWER_SHADERS})
    get_filename_component(SHADER_NAME ${SHADER} NAME)
    set(SPIRV_OUT ${VIEWER_SPIRV_DIR}/${SHADER_NAME}.spv)
    add_custom_command(
        OUTPUT ${SPIRV_OUT}
        COMMAND ${GLSLC} ${SHADER} -o ${SPIRV_OUT}
        DEPENDS ${SHADER}
        COMMENT "Compiling ${SHADER_NAME} -> SPIR-V"
    )
    list(APPEND VIEWER_SPIRV_OUTPUTS ${SPIRV_OUT})
endforeach()

add_custom_target(viewer_shaders DEPENDS ${VIEWER_SPIRV_OUTPUTS})

# --- Viewer executable ---

add_executable(viewer
    src/main.cpp
    src/vulkan_context.cpp
    src/renderer.cpp
    src/mesh.cpp
    src/camera.cpp
    src/ui.cpp
    src/fullerene_model.cpp
)

add_dependencies(viewer viewer_shaders)

target_include_directories(viewer PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_link_libraries(viewer PRIVATE
    Vulkan::Vulkan
    glfw
    glm::glm
    imgui
    GPUOpen::VulkanMemoryAllocator
    fullerenes
    ${LAPACK_LIBRARIES}
)

target_compile_definitions(viewer PRIVATE
    SHADER_DIR="${VIEWER_SPIRV_DIR}"
    GLM_FORCE_DEPTH_ZERO_TO_ONE
    VK_ENABLE_BETA_EXTENSIONS
)
