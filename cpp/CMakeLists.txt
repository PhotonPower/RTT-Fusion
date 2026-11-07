# C++ prototype of `rtt-fusion` (Fusion ADR draft point 9): port types and conversions, adapter
# interface; accepted against reference/*.result.json of the Python prototype.
cmake_minimum_required(VERSION 3.24)
project(rtt_fusion_cpp LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(nlohmann_json 3 CONFIG REQUIRED)

add_library(rtt_fusion_core STATIC src/ports.cpp src/diffraction_table.cpp)
target_include_directories(rtt_fusion_core PUBLIC include)
target_include_directories(rtt_fusion_core SYSTEM PRIVATE third_party/pocketfft)
if(MSVC)
  target_compile_options(rtt_fusion_core PRIVATE /W4)
else()
  target_compile_options(rtt_fusion_core PRIVATE -Wall -Wextra -Wpedantic)
endif()

enable_testing()
add_executable(test_reference tests/test_reference.cpp)
target_link_libraries(test_reference PRIVATE rtt_fusion_core nlohmann_json::nlohmann_json)
target_link_libraries(rtt_fusion_core PUBLIC nlohmann_json::nlohmann_json)
add_test(NAME reference COMMAND test_reference ${CMAKE_CURRENT_SOURCE_DIR}/../reference)
add_executable(test_diffraction_table tests/test_diffraction_table.cpp)
target_link_libraries(test_diffraction_table PRIVATE rtt_fusion_core)
add_test(NAME diffraction_table COMMAND test_diffraction_table ${CMAKE_CURRENT_SOURCE_DIR}/../reference)

# BWPM adapter (way 2: external C++ adapter against the public headers) – BWPM's C++ core from
# its repository at the pinned commit of pyproject.toml.
option(RTT_FUSION_BWPM_ADAPTER "Build the C++ BWPM adapter (fetches BWPM)" ON)
if(RTT_FUSION_BWPM_ADAPTER)
  include(FetchContent)
  FetchContent_Declare(bwpm GIT_REPOSITORY https://github.com/PhotonPower/BWPM.git
                       GIT_TAG ef2de941d226d364db000a977963b0eba24aea44)
  FetchContent_MakeAvailable(bwpm)
  add_library(rtt_fusion_bwpm_adapter STATIC adapters/bwpm/bwpm_adapter.cpp)
  target_include_directories(rtt_fusion_bwpm_adapter PUBLIC adapters/bwpm)
  target_link_libraries(rtt_fusion_bwpm_adapter PUBLIC rtt_fusion_core bwpm_core nlohmann_json::nlohmann_json)
  add_executable(test_bwpm_adapter tests/test_bwpm_adapter.cpp)
  target_link_libraries(test_bwpm_adapter PRIVATE rtt_fusion_bwpm_adapter)
  add_test(NAME bwpm_adapter COMMAND test_bwpm_adapter ${bwpm_SOURCE_DIR}/reference)
endif()
