#Make file
# please run ./configure before make
include Makefile.in
CC = $(CCINST) -fopenmp
LDFLAGS = -O3
OBJFLAGS = -c -Wall
DEBUGFLAGS = -g -Wall
PRECOMPILE =  -E -P
INC = -I include
PRE = precompiled
OBJ = obj
BIN = bin
DEBUG= debug
CFLAGS += 
GCCVERSION = $(shell gcc --version | grep ^gcc | sed 's/^.* //g')
DYNAMICFL= -D  PARFOR_DYNAMIC -D PAR_CHUNKSIZE=$(CHUNKSIZE)
GUIDEDFL= -D PARFOR_GUIDED
TASKFL= -D TASKLOOP_DEFINED

# (PARFOR_STATIC PARFOR_GUIDED PARFOR_DYNAMIC TASKLOOP)

SRCS=$(wildcard src/*.c)

OBJS=$(SRCS:.c=.o)

PROGS = $(patsubst %.c,%,$(SRCS))

all: bin

bin: $(OBJ)/%.o
	$(CC) $(LDFLAGS) -o $(basename $(notdir $<)) $<
$(OBJ)/%.o: $(SRCS)	
	$(CC) $(OBJFLAGS)  -D $(CAPABILITY) $(ONLINECORESFLAG) $(INC) -o $(OBJ)/$(basename $(notdir $<))_static.o   $<
	$(CC) $(OBJFLAGS)  -D $(CAPABILITY) $(ONLINECORESFLAG) $(INC) $(DYNAMICFL) -o $(OBJ)/$(basename $(notdir $<))_dynamic.o   $<
	$(CC) $(OBJFLAGS)  -D $(CAPABILITY) $(ONLINECORESFLAG) $(INC) $(GUIDEDFL) -o $(OBJ)/$(basename $(notdir $<))_guided.o   $<
	ifdef TASKLOOP_DEFINED
	$(CC) $(OBJFLAGS)   -D $(CAPABILITY) $(ONLINECORESFLAG) $(INC) $(TASKFL) -o $(OBJ)/$(basename $(notdir $<))_task.o   $<

debug: 
	$(CC) $(DEBUGFLAGS) $(INC) -D $(CAPABILITY) $(ONLINECORESFLAG) -o $(DEBUG)/$(basename $(notdir $<))_static $<
	$(CC) $(DEBUGFLAGS) $(INC) -D $(CAPABILITY) $(ONLINECORESFLAG) $(DYNAMICFL) -o $(DEBUG)/$(basename $(notdir $<))_dynamic $<
	$(CC) $(DEBUGFLAGS) $(INC) -D $(CAPABILITY) $(ONLINECORESFLAG) $(GUIDEDFL) -o $(DEBUG)/$(basename $(notdir $<))_guided $<
	ifdef TASKLOOP_DEFINED
	$(CC) $(DEBUGFLAGS) $(INC) -D $(CAPABILITY) $(ONLINECORESFLAG) $(TASKFL) -o $(DEBUG)/$(basename $(notdir $<))_task $<

fromprecompiled: precompile
	$(CC) $(LDFLAG) -o $(BIN)/$(basename $(notdir $<)) $<

precompile: src/%.c
	$(CC) $(PRECOMPILE) $(INC) -D $(CAPABILITY) $(ONLINECORESFLAG) -o  $(PRE)/$(basename $(notdir $<))_static.i  $<
	$(CC) $(PRECOMPILE) $(INC) -D $(CAPABILITY) $(ONLINECORESFLAG) $(DYNAMICFL) -o  $(PRE)/$(basename $(notdir $<))_dynamic.i  $<
	$(CC) $(PRECOMPILE) $(INC) -D $(CAPABILITY) $(ONLINECORESFLAG) $(GUIDEDFL) -o  $(PRE)/$(basename $(notdir $<))_guided.i  $<
	ifdef TASKLOOP_DEFINED
	$(CC) $(PRECOMPILE) $(INC) -D $(CAPABILITY) $(ONLINECORESFLAG) $(TASKFL) -o  $(PRE)/$(basename $(notdir $<))_task.i  $<

clean:
	rm -f $(OBJ)/*
	rm -f $(BIN)/*
	rm -f $(PRE)/*


