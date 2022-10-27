include ../include.mk

EXECUTABLE = ../../bin/optimizeComplexity

SOURCES    = optimizeComplexity.cpp \
			 ../../lib/Quaternion/POLYNOMIAL_4D.cpp \
			 ../../lib/Quaternion/QUATERNION.cpp \
			 ../../src/triangle.cpp \

OBJECTS = $(SOURCES:.cpp=.o)

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(WARNING) $(CXXFLAGS) $^ -o $@

.cpp.o:
	$(CXX) $(WARNING) $(CXXFLAGS) -MMD -MP -c $< -o $@

clean:
	rm -f *.o
