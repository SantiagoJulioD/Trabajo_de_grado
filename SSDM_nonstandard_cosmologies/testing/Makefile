# Compilador
CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Iinclude

# Directorios
SRC_DIR = src
OBJ_DIR = obj
INCLUDE_DIR = include
TEST_DIR = tests
BIN_DIR = bin

# Archivos fuente
SRC_SRCS = $(wildcard $(SRC_DIR)/*.cpp)
MAIN_SRC = main.cpp
TEST_SRCS = $(wildcard $(TEST_DIR)/*.cpp)

# Objetos
SRC_OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_SRCS))
MAIN_OBJ = $(OBJ_DIR)/main.o

# Ejecutables
TARGET = $(BIN_DIR)/main.out
TEST_BINS = $(patsubst $(TEST_DIR)/%.cpp,$(BIN_DIR)/%,$(TEST_SRCS))

# Regla principal
all: $(TARGET) $(TEST_BINS)

# Ejecutable principal (src + main.cpp)
$(TARGET): $(SRC_OBJS) $(MAIN_OBJ)
	mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Ejecutables de test (src + testX.cpp, pero sin main.o)
$(BIN_DIR)/%: $(TEST_DIR)/%.cpp $(SRC_OBJS)
	mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Compilación de objetos
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJ_DIR)/main.o: $(MAIN_SRC)
	mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Limpieza
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)


