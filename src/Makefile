# Compiler and linker configurations
CXX = g++
CXXFLAGS = $(shell root-config --cflags) -fPIC
LDFLAGS = $(shell root-config --libs)

# Target library
TARGET = libFileHandling.so
SOURCES = file_handling.cpp
OBJECTS = $(SOURCES:.cpp=.o)

# Default target
all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) -shared -o $@ $^ $(LDFLAGS)

# To obtain object files
%.o: %.cpp
	$(CXX) -c $< $(CXXFLAGS)

# Clean up
clean:
	rm -f $(OBJECTS) $(TARGET)

# Phony targets
.PHONY: all clean
