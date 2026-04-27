# For storing executables
OBJ_DIR = objs
# source files
SOURCES = $(wildcard src/*.cpp)
# object file path and name
EXECUTABLE = $(OBJ_DIR)/PFB_app

# command to compile the code
build:
	@echo "Compiling the code..."
	g++-15 -std=c++17 -Iinclude $(SOURCES) -o $(EXECUTABLE)

# command to run the executable
run: build
	@echo "Running the executable..."
	./$(EXECUTABLE)

# command to clean the build files
clean:
	@echo "Cleaning build files..."
	rm -rf $(OBJ_DIR)/*
