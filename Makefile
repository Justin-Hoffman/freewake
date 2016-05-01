#Compiler settings 
CC = gcc 
CPP = g++

#CPP = clang++-3.5 
#CPP = clang++-3.5

#CC = icc
#CPP = icpc
# Compile-time flags 
#CFLAGS = -Wall -Weffc++ -fPIC -std=c++11 -mtune=core-avx2 -msse4.2 -mavx2 -mpclmul -march=native -O3 -fprofile-use
CFLAGS = -Wall -Wextra -Weffc++ -fPIC -std=c++11 -g -O0 -fopenmp 
CPPFLAGS = -std=c++11 -L./ -lgcov -llapack -llapacke -lm -lgomp
#CPPFLAGS = -std=c++11 -L ./ -lpthread
INC = -I ./include -I ./src 
LDFLAGS =  -shared
SRCDIR  = ./src/

#gtest specific setup
GTEST_SRC_DIR = ./test/src/
GTEST_INC_DIR = ./test/inc/ -I /usr/src/gmock/gtest/include -I /usr/src/gmock/include
GTEST_CPP_FLAGS = $(CPPFLAGS) -lpthread -lgmock

# Where to find gtest_main.cc.
GTEST_MAIN_CC = ../../src/gtest_main.cc

# get all objects
CPP_FILES := $(wildcard ./src/*.cpp)
OBJ_FILES := $(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o)))

CPP_SRC_FILES = \
./src/horseshoelattice.cpp\
./src/liftingsurface.cpp\
./src/simulationmanager.cpp\
./src/tipfilament.cpp\
./src/vec3d.cpp\
./src/vortexlattice.cpp\
./src/vortexmath.cpp\


CPP_OBJ_FILES := $(addprefix obj/, $(notdir $(CPP_SRC_FILES:.cpp=.o)))


TEST_SRC_FILES = $(GTEST_SRC_DIR)gtest_main.cpp \
$(GTEST_SRC_DIR)horseshoelatticetest.cpp \
$(GTEST_SRC_DIR)liftingsurfacetest.cpp \
$(GTEST_SRC_DIR)simulationmanagertest.cpp \
$(GTEST_SRC_DIR)vec3dtest.cpp \
$(GTEST_SRC_DIR)vortexlatticetest.cpp \
$(GTEST_SRC_DIR)vortexmathtest.cpp \

#### TARGETS ####
all: freewake allTests check libfreewake.so

freewake: $(CPP_OBJ_FILES) obj/main.o
	$(CPP) \
$(CPP_OBJ_FILES) \
obj/main.o \
 $(CPPFLAGS) $(CFLAGS) $(INC) -o freewake

allTests: 
	$(CPP) \
$(TEST_SRC_FILES) \
$(CPP_OBJ_FILES) \
$(CPPFLAGS) $(GTEST_CPP_FLAGS$) $(INC) -I $(GTEST_INC_DIR) /usr/lib/libgtest.a -lpthread -o ./test/allTests

obj/%.o: src/%.cpp
	$(CC) $(CFLAGS) $(CPPFLAGS) $(INC) -c -o $@ $<

libfreewake.so: obj/liftingsurface.o obj/horseshoelattice.o obj/simulationmanager.o obj/vec3d.o obj/vortexlattice.o obj/vortexmath.o
	$(CPP) $(CPP_OBJ_FILES) $(CPPFLAGS) $(INC) $(LDFLAGS) -o ./libs/libfreewake.so 

check: 
	./test/allTests

clean:
	rm -rf ./src/*o freewake ./lib/libfreewake.so ./obj/*.o ./test/allTests
