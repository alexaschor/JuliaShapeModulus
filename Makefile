all: fast

CXX=g++-9
CXXFLAGS=-Wall -MMD -std=c++11 -fopenmp -I./lib/

# SOURCES_SRC := $(wildcard src/*.cpp)
SOURCES_SRC := src/main.cpp
SOURCES_LIB := $(wildcard lib/*/*.cpp)
SOURCES := $(SOURCES_SRC) $(SOURCES_LIB)

OBJECTS := $(patsubst %.cpp,build/%.o, $(SOURCES))
DEPENDS := $(patsubst build/%.o,build/%.d,$(OBJECTS))

# Warning level
WARNING := -Wall -Wextra

.PHONY: all clean

slow: CXXFLAGS += -g
fast: CXXFLAGS += -O3 -g

slow: run
fast: run

clean:
	rm -rf build run

# Linking the executable from the object files
run: $(OBJECTS)
	$(CXX) $(WARNING) $(CXXFLAGS) $^ -o $@

-include $(DEPENDS)

DIRGUARD := @mkdir -p $(dir $(DEPENDS) $(OBJECTS))

build/%.o: %.cpp Makefile
	$(DIRGUARD)
	$(CXX) $(WARNING) $(CXXFLAGS) -MMD -MP -c $< -o $@


