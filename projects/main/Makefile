include ../include.mk

EXECUTABLE = ../../bin/run

SOURCES    = main.cpp \
			 ../../lib/Quaternion/POLYNOMIAL_4D.cpp \
			 ../../lib/Quaternion/QUATERNION.cpp \

OBJECTS = $(SOURCES:.cpp=.o)

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(WARNING) $(CXXFLAGS) $^ -o $@

.cpp.o:
	$(CXX) $(WARNING) $(CXXFLAGS) -MMD -MP -c $< -o $@

clean:
	rm -f *.o
