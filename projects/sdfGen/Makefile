include ../include.mk

EXECUTABLE = ../../bin/sdfGen

SOURCES    = main.cpp \
			 ./makelevelset3.cpp

OBJECTS = $(SOURCES:.cpp=.o)

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(WARNING) $(CXXFLAGS) $^ -o $@

.cpp.o:
	$(CXX) $(WARNING) $(CXXFLAGS) -MMD -MP -c $< -o $@

clean:
	rm -f *.o
