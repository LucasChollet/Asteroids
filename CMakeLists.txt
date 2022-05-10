cmake_minimum_required(VERSION 3.12)
project(Asteroids
        VERSION 0.0.1
        LANGUAGES CXX)

include_directories(${PROJECT_SOURCE_DIR}/src)

option(SANITIZE "Compile with sanitizers" OFF)

set(CompilerWarnings -pedantic -Wall -Wextra -Wcast-align -Wcast-qual -Wdisabled-optimization
        -Wformat=2 -Winit-self -Wlogical-op -Wmissing-declarations -Wmissing-include-dirs -Wnoexcept
        -Woverloaded-virtual -Wredundant-decls -Wsign-promo -Wstrict-null-sentinel -Wstrict-overflow
        -Wswitch-default -Wundef -Weffc++ -Wconversion -Wrange-loop-construct)


set(Sanitizer -fsanitize=address -fsanitize=leak -fsanitize=undefined -fsanitize-address-use-after-scope)

add_executable(Asteroids src/main.cpp)

if (SANITIZE)
    message("Using sanitizers")
    target_compile_options(Asteroids PUBLIC ${CompilerWarnings} ${CompilerOptions} ${Sanitizer})
    target_link_options(Asteroids PUBLIC ${Sanitizer})
else()
    target_compile_options(Asteroids PUBLIC ${CompilerWarnings} ${CompilerOptions} -Os)
endif (SANITIZE)

find_package(glm REQUIRED)
include_directories(${glm_INCLUDE_DIRS})
link_directories(${glm_LIBRARY_DIRS})

find_package(gcem REQUIRED)
include_directories(${gcem_INCLUDE_DIRS})
link_directories(${gcem_LIBRARY_DIRS})

set_property(TARGET Asteroids PROPERTY CXX_STANDARD 20)
set_property(TARGET Asteroids PROPERTY CMAKE_CXX_STANDARD_REQUIRED ON)

if (BUILD_DOC)
    find_package(Doxygen REQUIRED dot)
    include(CMakeDoxygen.txt)

    doxygen_add_docs(Doxygen ALL
            src README.md
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            COMMENT "Generating API documentation with Doxygen"
            )
endif (BUILD_DOC)
