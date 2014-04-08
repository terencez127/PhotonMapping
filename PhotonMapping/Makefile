#
# Generic makefile for Mac OSX an example program using the cprocessing library
#
OSTYPE := $(shell uname)
CFLAGS = -g -fPIC

LDFLAGS = -lGL -lglut -lGLU
ifeq ("$(OSTYPE)", "Darwin")
    LDFLAGS = -framework Carbon -framework OpenGL -framework GLUT
    export MACOSX_DEPLOYMENT_TARGET=10.8
endif

CPROCESSINGPATH=/Users/terence/dev/cprocessing

CPROCESSINGLIB=$(CPROCESSINGPATH)/objects

SOURCES= $(wildcard *.cpp)

OBJECTS=$(patsubst %.cpp, %.o, $(SOURCES))

TARGET=$(lastword $(subst /, ,$(realpath .)))

all: $(TARGET)
	@echo $(TARGET) built

$(TARGET): $(OBJECTS) $(CPROCESSINGLIB)/libcprocessing.a
	@clang++ -std=c++11 -stdlib=libc++ -Werror -mmacosx-version-min=10.8 -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.8.sdk $(CFLAGS) -o $@ $(OBJECTS) -L$(CPROCESSINGLIB) $(LDFLAGS) -lcprocessing

clean:
	rm -rf $(OBJECTS) $(TARGET)

%.o: %.cpp $(CPROCESSINGPATH)/cprocessing.hpp
	@echo compiling $<
	@clang++ -std=c++11 -stdlib=libc++ -Werror $(CFLAGS) -c $< -o $@ -I$(CPROCESSINGPATH)
