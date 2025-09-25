# Compilador
CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Iinclude

# Archivos fuente y objeto
SRC_DIR = src
OBJ_DIR = obj
INCLUDE_DIR = include

SRCS = $(wildcard $(SRC_DIR)/*.cpp) main.cpp
OBJS = $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(notdir $(SRCS)))

# Nombre del ejecutable
TARGET = main.out

# Regla principal
all: $(TARGET)

# Compilación del ejecutable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Compilar archivos .cpp a .o
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJ_DIR)/main.o: main.cpp
	mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Limpieza
clean:
	rm -rf $(OBJ_DIR) *.out
