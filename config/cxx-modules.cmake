# Default to OFF
set(CADO_USE_STD_MODULE)

# Check for CMake 3.30+ and Ninja generator
if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.30" AND CMAKE_GENERATOR MATCHES "Ninja")
    message(STATUS "Checking support for C++23 'import std;'...")
    
    # Create a minimal test project to verify the toolchain
    file(WRITE "${CMAKE_BINARY_DIR}/module_test/CMakeLists.txt" [=[
        cmake_minimum_required(VERSION 3.30)
        project(module_test CXX)
        set(CMAKE_CXX_STANDARD 23)
        set(CMAKE_CXX_STANDARD_REQUIRED ON)
        set(CMAKE_CXX_EXTENSIONS ON)
        add_executable(test_std main.cpp)
        set_target_properties(test_std PROPERTIES 
            CXX_STANDARD 23 
            CXX_MODULE_STD ON
        )
    ]=])
    
    file(WRITE "${CMAKE_BINARY_DIR}/module_test/main.cpp" [=[
        import std;
        int main() { std::cout << "Modules work!\n"; return 0; }
    ]=])

    try_compile(HAVE_CXX_MODULE_STD
        PROJECT module_test
        SOURCE_DIR "${CMAKE_BINARY_DIR}/module_test"
    )

    if (HAVE_CXX_MODULE_STD)
        set(CADO_USE_STD_MODULE 1)
        message(STATUS "Checking support for C++23 'import std;' -- Success")
    else()
        message(STATUS "Checking support for C++23 'import std;' -- Failed (Compiler/Stdlib mismatch)")
    endif()
endif()

function(cado_enable_std_module target)
    if(CADO_USE_STD_MODULE)
        set_target_properties(${target} PROPERTIES
            CXX_STANDARD 23
            CXX_MODULE_STD ON
        )
    endif()
endfunction()
