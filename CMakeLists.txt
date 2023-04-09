cmake_minimum_required(VERSION 3.20)
project (PairwiseSequenceAlignment)

set(CMAKE_CXX_STANDARD 17)

add_executable(PairwiseSequenceAlignment DynamicProgramming.cpp)
add_executable(RecursiveAlignment RecursiveFunction.cpp)
