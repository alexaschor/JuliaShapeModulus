include ../include.mk

EXECUTABLE = ../../bin/rootMorph

SOURCES    = rootMorph.cpp \
			 ../../lib/Quaternion/QUATERNION.cpp \
			 ../../lib/Quaternion/POLYNOMIAL_4D.cpp \
			 ../../src/triangle.cpp

OBJECTS = $(SOURCES:.cpp=.o)

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(WARNING) $(CXXFLAGS) $^ -o $@

.cpp.o:
	$(CXX) $(WARNING) $(CXXFLAGS) -MMD -MP -c $< -o $@

clean:
	rm -f *.o
