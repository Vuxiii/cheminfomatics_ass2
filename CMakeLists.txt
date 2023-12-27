cmake_minimum_required(VERSION 3.10)

project(CMakeTestProject CXX)

option(WITH_PYTHON_EXT "Whether to build it as a Python extension as well or not." ON)

if(NOT CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

if(CMAKE_COMPILER_IS_GNUCXX OR ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang"))
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
endif()

find_package(mod REQUIRED)

if(WITH_PYTHON_EXT)
    # Boost.Python
    # -------------------------------------------------------------------------
    set(v 1.64.0)
    foreach(PY 3 37 38 39 310 311 312 313 314 315 316 317 318 319)
        set(lib "python${PY}")
        find_package(Boost ${v} QUIET COMPONENTS ${lib})
        if(Boost_FOUND)
            find_package(Boost ${v} COMPONENTS ${lib})
            set(PYTHON_TARGET ${lib})
            break()
        endif()
    endforeach()
    if(NOT Boost_FOUND)
        find_package(Boost ${v} REQUIRED COMPONENTS python3)
        message(FATAL_ERROR "Could not find Boost.Python for Python 3. Tried 'python' wih suffixes 3, 37, 38, 39, and 310 to 319.")
    endif()

    find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
endif()

# add automatically determined dependencies outside the build dir to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# As executable
add_executable(doStuff main.cpp)
target_link_libraries(doStuff mod::libmod)
install(TARGETS doStuff RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX})

# As a Python extension
if(WITH_PYTHON_EXT)
    add_library(pydoStuff MODULE main.cpp)
    set_target_properties(pydoStuff PROPERTIES PREFIX "") # so it doesn't get the "lib" prefix
    target_compile_definitions(pydoStuff PRIVATE AS_PYTHON_EXTENSION)
    target_link_libraries(pydoStuff mod::libmod Boost::${PYTHON_TARGET} Python3::Python)
    install(TARGETS pydoStuff LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX})
endif()
