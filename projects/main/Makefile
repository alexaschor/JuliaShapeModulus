include ../include.mk

EXECUTABLE = ../../bin/run

SOURCES    = main.cpp \
			 ../../lib/Quaternion/POLYNOMIAL_4D.cpp \
			 ../../lib/Quaternion/QUATERNION.cpp \

OBJECTS = $(SOURCES:.cpp=.o)

all: $(SOURCES) $(EXECUTABLE) script-link

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(WARNING) $(CXXFLAGS) $^ -o $@

.cpp.o:
	$(CXX) $(WARNING) $(CXXFLAGS) -MMD -MP -c $< -o $@

script-link:
	@ln -fs ../projects/main/prun.py ../../bin/prun

clean:
	rm -f *.o
