#Make file
# please run ./configure before make
include Makefile.in
CC =$(CCINST) -fopenmp
LDFLAGS = -O3 
EM= -lenergymodule
OBJFLAGS = -c -Wall
DEBUGFLAGS = -g -Wall
PRECOMPILE =  -E -P
INC = -I ./include/ -I../Energymonitorlibrary/include/
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
TASKFL= -D TASKLOOP_DEFINED -D NUM_TASKS=$(NUMTASKS)
ONLINECORESFLAG=-D ONLINECORES=$(ONLINE_CORES)


# (PARFOR_STATIC PARFOR_GUIDED PARFOR_DYNAMIC TASKLOOP)

SRCS=$(wildcard src/*.c)
DEBUGS=$(wildcard src/*.c)

OBJS=$(SRCS:.c=.o)

PROGS = $(patsubst %.c,%,$(SRCS))

.phony: $(SRCS) 

all: $(SRCS) 
debug: $(DEBUGS)

# ifeq ($(TASKLOOP_DEFINED), yes)
# bin: $(SRCS)  taskobj
# else
# bin: $(SRCS)
# endif
# 	$(CC) -D $(CAPABILITY) $(ONLINECORESFLAG) $(INC)  $(LDFLAGS) $<  $(EM) -o $(BIN)/$(basename $(notdir $<))_static  
# 	$(CC) -D $(CAPABILITY) $(ONLINECORESFLAG) $(INC) $(LDFLAGS) $(DYNAMICFL) $< $(EM) -o $(BIN)/$(basename $(notdir $<))_dynamic   
# 	$(CC) -D $(CAPABILITY) $(ONLINECORESFLAG) $(INC) $(LDFLAGS) $(GUIDEDFL) $< $(EM) -o $(BIN)/$(basename $(notdir $<))_guided  
# taskobj: $(SRCS) 
# 	$(CC) -D $(CAPABILITY) $(ONLINECORESFLAG) $(INC) $(LDFLAGS) $(TASKFL) $< $(EM) -o $(OBJ)/$(basename $(notdir $<))_task


#######	$(CC) $(LDFLAGS) -o $(basename $(notdir $<)) $
ifeq ($(TASKLOOP_DEFINED), yes)
$(SRCS):
	$(CC) -D $(CAPABILITY) $(ONLINECORESFLAG) $(INC) $(LDFLAGS) $(TASKFL) $@ $(EM) -o $(BIN)/$(basename $(notdir $@))_task
else
$(SRCS): 
endif
	$(CC) -D $(CAPABILITY) $(ONLINECORESFLAG) $(INC)  $(LDFLAGS) $@  $(EM) -o $(BIN)/$(basename $(notdir $@))_static  
	$(CC) -D $(CAPABILITY) $(ONLINECORESFLAG) $(INC) $(LDFLAGS) $(DYNAMICFL) $@ $(EM) -o $(BIN)/$(basename $(notdir $@))_dynamic 
	$(CC) -D $(CAPABILITY) $(ONLINECORESFLAG) $(INC) $(LDFLAGS) $(GUIDEDFL) $@ $(EM) -o $(BIN)/$(basename $(notdir $@))_guided  


dynamicwithchunkSize:
	$(CC) -D $(CAPABILITY) $(ONLINECORESFLAG) $(INC) $(LDFLAGS) $(DYNAMICFL) $(TARGET) $(EM) -o $(BIN)/$(basename $(notdir $(TARGET)))_dynamic_$(CHUNKSIZE) 


# ifeq ($(TASKLOOP_DEFINED), yes )
# $(DEBUGS):
# 	$(CC) $(DEBUGFLAGS) $(INC) -D $(CAPABILITY) $(ONLINECORESFLAG) $(TASKFL)  $(LDFLAGS) $@ $(EM) -o $(DEBUG)/$(basename $(notdir $@))_task
# else
# $(DEBUGS):
# endif
# 	$(CC) $(DEBUGFLAGS) $(INC) -D $(CAPABILITY) $(ONLINECORESFLAG)  $(LDFLAGS) $@ $(EM) -o $(DEBUG)/$(basename $(notdir $@))_static
# 	$(CC) $(DEBUGFLAGS) $(INC) -D $(CAPABILITY) $(ONLINECORESFLAG) $(DYNAMICFL) $(LDFLAGS) $@ $(EM) -o $(DEBUG)/$(basename $(notdir $<))_dynamic
# 	$(CC) $(DEBUGFLAGS) $(INC) -D $(CAPABILITY) $(ONLINECORESFLAG) $(GUIDEDFL)  $(LDFLAGS) $@ $(EM) -o $(DEBUG)/$(basename $(notdir $<))_guided


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


