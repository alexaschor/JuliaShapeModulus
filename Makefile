CXX=g++
CXXFLAGS=-g -Wall -MMD -std=c++11

SOURCES := $(wildcard *.cpp)
OBJECTS := $(patsubst %.cpp,%.o,$(SOURCES))
DEPENDS := $(patsubst %.cpp,%.d,$(SOURCES))

# Warning level
WARNING := -Wall -Wextra

.PHONY: all clean

all: run

clean:
	$(RM) $(OBJECTS) $(DEPENDS) run

# Linking the executable from the object files
run: $(OBJECTS)
	$(CXX) $(WARNING) $(CXXFLAGS) $^ -o $@

-include $(DEPENDS)

%.o: %.cpp Makefile
	$(CXX) $(WARNING) $(CXXFLAGS) -MMD -MP -c $< -o $@
