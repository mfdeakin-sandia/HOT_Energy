# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.6

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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
CMAKE_COMMAND = /usr/bin/cmake3

# The command to remove a file.
RM = /usr/bin/cmake3 -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /ascldap/users/scmousl/HOT_Energy

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /ascldap/users/scmousl/HOT_Energy/build

# Include any dependencies generated for this target.
include CMakeFiles/wassertest.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/wassertest.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/wassertest.dir/flags.make

CMakeFiles/wassertest.dir/src/wasserstein/wasserstein.cpp.o: CMakeFiles/wassertest.dir/flags.make
CMakeFiles/wassertest.dir/src/wasserstein/wasserstein.cpp.o: ../src/wasserstein/wasserstein.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/ascldap/users/scmousl/HOT_Energy/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/wassertest.dir/src/wasserstein/wasserstein.cpp.o"
	/usr/bin/clang++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/wassertest.dir/src/wasserstein/wasserstein.cpp.o -c /ascldap/users/scmousl/HOT_Energy/src/wasserstein/wasserstein.cpp

CMakeFiles/wassertest.dir/src/wasserstein/wasserstein.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/wassertest.dir/src/wasserstein/wasserstein.cpp.i"
	/usr/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /ascldap/users/scmousl/HOT_Energy/src/wasserstein/wasserstein.cpp > CMakeFiles/wassertest.dir/src/wasserstein/wasserstein.cpp.i

CMakeFiles/wassertest.dir/src/wasserstein/wasserstein.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/wassertest.dir/src/wasserstein/wasserstein.cpp.s"
	/usr/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /ascldap/users/scmousl/HOT_Energy/src/wasserstein/wasserstein.cpp -o CMakeFiles/wassertest.dir/src/wasserstein/wasserstein.cpp.s

CMakeFiles/wassertest.dir/src/wasserstein/wasserstein.cpp.o.requires:

.PHONY : CMakeFiles/wassertest.dir/src/wasserstein/wasserstein.cpp.o.requires

CMakeFiles/wassertest.dir/src/wasserstein/wasserstein.cpp.o.provides: CMakeFiles/wassertest.dir/src/wasserstein/wasserstein.cpp.o.requires
	$(MAKE) -f CMakeFiles/wassertest.dir/build.make CMakeFiles/wassertest.dir/src/wasserstein/wasserstein.cpp.o.provides.build
.PHONY : CMakeFiles/wassertest.dir/src/wasserstein/wasserstein.cpp.o.provides

CMakeFiles/wassertest.dir/src/wasserstein/wasserstein.cpp.o.provides.build: CMakeFiles/wassertest.dir/src/wasserstein/wasserstein.cpp.o


# Object files for target wassertest
wassertest_OBJECTS = \
"CMakeFiles/wassertest.dir/src/wasserstein/wasserstein.cpp.o"

# External object files for target wassertest
wassertest_EXTERNAL_OBJECTS =

wassertest: CMakeFiles/wassertest.dir/src/wasserstein/wasserstein.cpp.o
wassertest: CMakeFiles/wassertest.dir/build.make
wassertest: /usr/lib64/libmpfr.so
wassertest: /usr/lib64/libgmp.so
wassertest: /usr/local/lib64/libCGAL_Core.so.13.0.1
wassertest: /usr/local/lib64/libCGAL.so.13.0.1
wassertest: /usr/lib64/libboost_thread-mt.so
wassertest: /usr/lib64/libboost_system-mt.so
wassertest: /usr/lib64/libboost_chrono-mt.so
wassertest: /usr/lib64/libboost_date_time-mt.so
wassertest: /usr/lib64/libboost_atomic-mt.so
wassertest: /usr/lib64/libboost_thread-mt.so
wassertest: /usr/lib64/libboost_system-mt.so
wassertest: /usr/lib64/libboost_chrono-mt.so
wassertest: /usr/lib64/libboost_date_time-mt.so
wassertest: /usr/lib64/libboost_atomic-mt.so
wassertest: /usr/lib64/libboost_thread-mt.so
wassertest: /usr/lib64/libboost_system-mt.so
wassertest: /usr/lib64/libboost_chrono-mt.so
wassertest: /usr/lib64/libboost_date_time-mt.so
wassertest: /usr/lib64/libboost_atomic-mt.so
wassertest: libwasserstein.a
wassertest: /usr/local/lib64/libCGAL_Core.so.13.0.1
wassertest: /usr/local/lib64/libCGAL.so.13.0.1
wassertest: /usr/lib64/libboost_thread-mt.so
wassertest: /usr/lib64/libboost_system-mt.so
wassertest: /usr/lib64/libboost_chrono-mt.so
wassertest: /usr/lib64/libboost_date_time-mt.so
wassertest: /usr/lib64/libboost_atomic-mt.so
wassertest: /usr/lib64/libboost_thread-mt.so
wassertest: /usr/lib64/libboost_system-mt.so
wassertest: /usr/lib64/libboost_chrono-mt.so
wassertest: /usr/lib64/libboost_date_time-mt.so
wassertest: /usr/lib64/libboost_atomic-mt.so
wassertest: /usr/lib64/libboost_thread-mt.so
wassertest: /usr/lib64/libboost_system-mt.so
wassertest: /usr/lib64/libboost_chrono-mt.so
wassertest: /usr/lib64/libboost_date_time-mt.so
wassertest: /usr/lib64/libboost_atomic-mt.so
wassertest: /usr/lib64/libmpfr.so
wassertest: /usr/lib64/libgmp.so
wassertest: /usr/local/lib64/libCGAL_Core.so.13.0.1
wassertest: /usr/local/lib64/libCGAL.so.13.0.1
wassertest: /usr/lib64/libboost_thread-mt.so
wassertest: /usr/lib64/libboost_system-mt.so
wassertest: /usr/lib64/libboost_chrono-mt.so
wassertest: /usr/lib64/libboost_date_time-mt.so
wassertest: /usr/lib64/libboost_atomic-mt.so
wassertest: /usr/lib64/libboost_thread-mt.so
wassertest: /usr/lib64/libboost_system-mt.so
wassertest: /usr/lib64/libboost_chrono-mt.so
wassertest: /usr/lib64/libboost_date_time-mt.so
wassertest: /usr/lib64/libboost_atomic-mt.so
wassertest: CMakeFiles/wassertest.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/ascldap/users/scmousl/HOT_Energy/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable wassertest"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/wassertest.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/wassertest.dir/build: wassertest

.PHONY : CMakeFiles/wassertest.dir/build

CMakeFiles/wassertest.dir/requires: CMakeFiles/wassertest.dir/src/wasserstein/wasserstein.cpp.o.requires

.PHONY : CMakeFiles/wassertest.dir/requires

CMakeFiles/wassertest.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/wassertest.dir/cmake_clean.cmake
.PHONY : CMakeFiles/wassertest.dir/clean

CMakeFiles/wassertest.dir/depend:
	cd /ascldap/users/scmousl/HOT_Energy/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /ascldap/users/scmousl/HOT_Energy /ascldap/users/scmousl/HOT_Energy /ascldap/users/scmousl/HOT_Energy/build /ascldap/users/scmousl/HOT_Energy/build /ascldap/users/scmousl/HOT_Energy/build/CMakeFiles/wassertest.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/wassertest.dir/depend
