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
CMAKE_SOURCE_DIR = /Users/aryan/gd/USI/hpc/project/code/build/v0

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/aryan/gd/USI/hpc/project/code/build/v0

# Utility rule file for sign.

# Include the progress variables for this target.
include CMakeFiles/sign.dir/progress.make

CMakeFiles/sign:
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/aryan/gd/USI/hpc/project/code/build/v0/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Digitally signing main"
	/opt/local/bin/cmake -E echo '' && /opt/local/bin/cmake -E echo '***************************************************************************' && /opt/local/bin/cmake -E echo '** Error: No Mac OSX developer certificate specified' && /opt/local/bin/cmake -E echo '** Please reconfigure with -DOSX_CERTIFICATE_NAME="<...>"' && /opt/local/bin/cmake -E echo '***************************************************************************' && /opt/local/bin/cmake -E echo ''

sign: CMakeFiles/sign
sign: CMakeFiles/sign.dir/build.make
.PHONY : sign

# Rule to build all files generated by this target.
CMakeFiles/sign.dir/build: sign
.PHONY : CMakeFiles/sign.dir/build

CMakeFiles/sign.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/sign.dir/cmake_clean.cmake
.PHONY : CMakeFiles/sign.dir/clean

CMakeFiles/sign.dir/depend:
	cd /Users/aryan/gd/USI/hpc/project/code/build/v0 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/aryan/gd/USI/hpc/project/code/build/v0 /Users/aryan/gd/USI/hpc/project/code/build/v0 /Users/aryan/gd/USI/hpc/project/code/build/v0 /Users/aryan/gd/USI/hpc/project/code/build/v0 /Users/aryan/gd/USI/hpc/project/code/build/v0/CMakeFiles/sign.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/sign.dir/depend

