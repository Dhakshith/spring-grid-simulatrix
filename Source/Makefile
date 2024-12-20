OBJS = main.cpp
CC = em++
COMPILER_FLAGS = -sUSE_SDL=2 -sUSE_SDL_TTF=2
LINKER_FLAGS = -lm --preload-file Resources -s EXPORTED_RUNTIME_METHODS=["ccall","stringToNewUTF8"] -s EXPORTED_FUNCTIONS=_main,_setJsonString
OBJ_NAME = ../docs/index.js

all : $(OBJS)
	$(CC) $(OBJS) $(COMPILER_FLAGS) $(LINKER_FLAGS) -o $(OBJ_NAME)
