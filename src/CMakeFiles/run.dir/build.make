# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.0

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/local/bin/cmake

# The command to remove a file.
RM = /opt/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/aryan/gd/Projects/Multi-Asset-American-Options-FEA/src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/aryan/gd/Projects/Multi-Asset-American-Options-FEA/src

# Utility rule file for run.

# Include the progress variables for this target.
include CMakeFiles/run.dir/progress.make

CMakeFiles/run: main
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/aryan/gd/Projects/Multi-Asset-American-Options-FEA/src/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Run main with Debug configuration"
	./main

run: CMakeFiles/run
run: CMakeFiles/run.dir/build.make
.PHONY : run

# Rule to build all files generated by this target.
CMakeFiles/run.dir/build: run
.PHONY : CMakeFiles/run.dir/build

CMakeFiles/run.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/run.dir/cmake_clean.cmake
.PHONY : CMakeFiles/run.dir/clean

CMakeFiles/run.dir/depend:
	cd /Users/aryan/gd/Projects/Multi-Asset-American-Options-FEA/src && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/aryan/gd/Projects/Multi-Asset-American-Options-FEA/src /Users/aryan/gd/Projects/Multi-Asset-American-Options-FEA/src /Users/aryan/gd/Projects/Multi-Asset-American-Options-FEA/src /Users/aryan/gd/Projects/Multi-Asset-American-Options-FEA/src /Users/aryan/gd/Projects/Multi-Asset-American-Options-FEA/src/CMakeFiles/run.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/run.dir/depend

