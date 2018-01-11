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
ifdef SET_CHUNKSIZE
	CHUNKSIZE=$(SET_CHUNKSIZE)
endif
DYNAMICFL= -D  PARFOR_DYNAMIC -D PAR_CHUNKSIZE=$(CHUNKSIZE)
GUIDEDFL= -D PARFOR_GUIDED -D PAR_CHUNKSIZE=$(CHUNKSIZE)
TASKFL= -D TASKLOOP_DEFINED


# (PARFOR_STATIC PARFOR_GUIDED PARFOR_DYNAMIC TASKLOOP)

SRCS=$(wildcard src/*.c)

OBJS=$(SRCS:.c=.o)

PROGS = $(patsubst %.c,%,$(SRCS))

all: bin

ifeq ($(TASKLOOP_DEFINED), yes)
bin: $(OBJ)/%.o  taskobj
endif
bin: $(OBJ)/%.o  
	$(CC) $(LDFLAGS) -o $(basename $(notdir $<)) $<

$(OBJ)/%.o: $(SRCS)	
	$(CC) $(OBJFLAGS)  -D $(CAPABILITY) $(ONLINECORESFLAG) $(INC) -o $(OBJ)/$(basename $(notdir $<))_static.o   $<
	$(CC) $(OBJFLAGS)  -D $(CAPABILITY) $(ONLINECORESFLAG) $(INC) $(DYNAMICFL) -o $(OBJ)/$(basename $(notdir $<))_dynamic.o   $<
	$(CC) $(OBJFLAGS)  -D $(CAPABILITY) $(ONLINECORESFLAG) $(INC) $(GUIDEDFL) -o $(OBJ)/$(basename $(notdir $<))_guided.o   $<
taskobj: $(SRCS)
	$(CC) $(OBJFLAGS)   -D $(CAPABILITY) $(ONLINECORESFLAG) $(INC) $(TASKFL) -o $(OBJ)/$(basename $(notdir $<))_task.o   $<

ifeq ($(TASKLOOP_DEFINED), yes )
debug:
	$(CC) $(DEBUGFLAGS) $(INC) -D $(CAPABILITY) $(ONLINECORESFLAG) $(TASKFL) -o $(DEBUG)/$(basename $(notdir $<))_task $<
else
debug:
endif
	$(CC) $(DEBUGFLAGS) $(INC) -D $(CAPABILITY) $(ONLINECORESFLAG) -o $(DEBUG)/$(basename $(notdir $<))_static $<
	$(CC) $(DEBUGFLAGS) $(INC) -D $(CAPABILITY) $(ONLINECORESFLAG) $(DYNAMICFL) -o $(DEBUG)/$(basename $(notdir $<))_dynamic $<
	$(CC) $(DEBUGFLAGS) $(INC) -D $(CAPABILITY) $(ONLINECORESFLAG) $(GUIDEDFL) -o $(DEBUG)/$(basename $(notdir $<))_guided $<


fromprecompiled: precompile
	$(CC) $(LDFLAG) -o $(BIN)/$(basename $(notdir $<)) $<

ifeq ( $(TASKLOOP_DEFINED), yes )
precompile: src/%.c
	$(CC) $(PRECOMPILE) $(INC) -D $(CAPABILITY) $(ONLINECORESFLAG) $(TASKFL) -o  $(PRE)/$(basename $(notdir $<))_task.i  $<
else
precompile: src/%.c
endif
	$(CC) $(PRECOMPILE) $(INC) -D $(CAPABILITY) $(ONLINECORESFLAG) -o  $(PRE)/$(basename $(notdir $<))_static.i  $<
	$(CC) $(PRECOMPILE) $(INC) -D $(CAPABILITY) $(ONLINECORESFLAG) $(DYNAMICFL) -o  $(PRE)/$(basename $(notdir $<))_dynamic.i  $<
	$(CC) $(PRECOMPILE) $(INC) -D $(CAPABILITY) $(ONLINECORESFLAG) $(GUIDEDFL) -o  $(PRE)/$(basename $(notdir $<))_guided.i  $<

clean:
	rm -f $(OBJ)/*
	rm -f $(BIN)/*
	rm -f $(PRE)/*
	rm -f $(DEBUG)/*


