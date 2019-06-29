cmake_minimum_required(VERSION 3.6)

#find_package (Eigen3 3.3 REQUIRED )
project(Swapodopolis)


if( WIN32 )
        set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
        set(BUILD_SHARED_LIBS TRUE)

        set(Boost_USE_STATIC_LIBS        ON) # only find static libs
        set(Boost_USE_MULTITHREADED      ON)
        set(Boost_USE_STATIC_RUNTIME    OFF)

        set(BOOST_INCLUDEDIR C:/work/modular-boost )
        set(BOOST_LIBRARYDIR C:/work/modular-boost/stage/lib )
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /D_CRT_SECURE_NO_WARNINGS")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MT")
else()
        add_definitions(-DBOOST_LOG_DYN_LINK)
        set(CMAKE_CXX_FLAGS "-std=c++14 -frtti -Wno-switch")

        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb3")
        #set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall  -pedantic")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-sign-compare")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-variable")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-missing-braces")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-misleading-indentation")

endif()



find_package(Boost REQUIRED COMPONENTS log system timer serialization)

include_directories(${Boost_INCLUDE_DIRS})
include_directories(/usr/include/eigen3/)

include_directories(.)
include_directories(deps/QuantLib )
include_directories(deps/CandyPretty/include)
include_directories(/home/dest/repo/eigen-git-mirror/)


function(swapodopolis_add_example exe)
        add_executable(${exe} ${ARGN})
        target_link_libraries(${exe} ${Boost_LIBRARIES})
        target_link_libraries(${exe} QuantLib.so )
        target_link_libraries(${exe} pthread )
endfunction()

swapodopolis_add_example(proc proc.cpp)

