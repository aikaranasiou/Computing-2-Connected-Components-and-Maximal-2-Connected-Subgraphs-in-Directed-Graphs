# ##################################################################
# 	Make file to Compute Two Vertex Strongly Connected Component   #	
####################################################################
#
TARGET  = 2V-CHILP
SOURCES = 2V-CHILP.cpp rfw_timer.cpp
#
GCC_NAME    = g++
GCC_FLAGS   = -g -Wall -O3 -ansi -pedantic -fpermissive -std=c++11
GCC_LIBS    = -lm -L/usr/lib/ -lrt  -static-libstdc++
GCC_DEFINES = -DRFW_RUSAGE
GCC_OBJECTS = $(SOURCES:.cpp=.o)

VCC_NAME    = cl 
VCC_FLAGS   = -O2 /W4
VCC_DEFINES = -DWIN32 -DNDEBUG -D_CONSOLE 
VCC_LIBS    = -lm
VCC_OBJECTS = $(SOURCES:.cpp=.obj)

CC_NAME     = CC
CC_FLAGS    = -O3 -OPT:Olimit_opt=on -LANG:std
#CC_FLAGS    = -LANG:std
CC_DEFINES  = -DRFW_RUSAGE
CC_LIBS     = -lm -L/usr/lib32/mips3/
CC_OBJECTS  = $(SOURCES:.cpp=.o)

#
# CHANGE THESE LINES TO USE YOUR FAVORITE COMPILER
CCC      = $(GCC_NAME)
FLAGS    = $(GCC_FLAGS)
LIBS     = $(GCC_LIBS)
DEFINES  = $(GCC_DEFINES)
OBJECTS  = $(GCC_OBJECTS)
INCLUDES = -I.

.SUFFIXES: .cpp

# $< is replaced with name of the input file
# $@ is replaced with the name out the output file.

.cpp.o:
	$(CCC) $(DEFINES) $(FLAGS) -c $< -o $@

.cpp.obj:
	$(CCC) $(DEFINES) $(FLAGS) -c $< -o $@

$(TARGET): $(OBJECTS)
	$(CCC) $(FLAGS) $(DEFINES) $(INCLUDES) $(OBJECTS) $(LIBS) -o $(TARGET)


all: make clean

clean: 
	rm -f $(OBJECTS) #$(TARGET) *~
