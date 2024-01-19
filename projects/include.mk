# Settings for all project Makefiles

CXX=g++
CXXFLAGS=-Wall -MMD -g -Ofast -std=c++17 -fopenmp -I../../lib -I../../src/ -I../../
LD_PATH_VAR=LD_LIBRARY_PATH
