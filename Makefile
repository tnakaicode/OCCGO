TARGETS = sdl2_sample

all: $(TARGETS)

SDL_PREFIX	= /mingw64
SDL_CONFIG	= $(SDL_PREFIX)/bin/sdl2-config
CG_LIBS		= 

CROSS_COMPILE	= /mingw64/bin/
CC		= $(CROSS_COMPILE)gcc
CXX		= $(CROSS_COMPILE)g++

CFLAGS 		= -g -Wall `/bin/sh $(SDL_CONFIG) --cflags`
CXXFLAGS	= -g -Wall `/bin/sh $(SDL_CONFIG) --cflags`
LDFLAGS		= `/bin/sh $(SDL_CONFIG) --libs`	-Wl,-rpath,$(SDL_PREFIX)/lib
LIBS		= -lopengl32 -lglu32 -lm

clean:
	rm -f *.o *.a *~ $(TARGETS)

sdl2_sample: sdl2_sample.o
	$(CXX) -o $@ $^ $(LDFLAGS) $(LIBS)