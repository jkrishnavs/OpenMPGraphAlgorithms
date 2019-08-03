#Make file
# please run ./configure before make
include Makefile.in
CC =$(CCINST) -fopenmp
CPP =$(CPPINST) -fopenmp
LDFLAGS = -O3 
EM= -lenergymodule
OBJFLAGS = -c -Wall
DEBUGFLAGS = -g -Wall
PRECOMPILE =  -E -P
INCWITOUTBIGLITTLE=-I ./include/
INC = -I ./include/ -I../Energymonitorlibrary/include/ -D BIGLITTLE_H
PRE = precompiled
OBJ = obj
BIN = bin
DEBUG= debug
CFLAGS +=
SPRNG_BIN=(convert.tmp pi-simple.tmp seed.tmp seed-simple.tmp simple-simple.tmp spawn.tmp sprng.tmp sprng-simple.tmp)
GCCVERSION = $(shell gcc --version | grep ^gcc | sed 's/^.* //g')
ifdef SET_CHUNKSIZE
	CHUNKSIZE=$(SET_CHUNKSIZE)
endif
DYNAMICFL= -D  PARFOR_DYNAMIC -D PAR_CHUNKSIZE=$(CHUNKSIZE)
GUIDEDFL= -D PARFOR_GUIDED -D PAR_CHUNKSIZE=$(CHUNKSIZE)
TASKFL= -D TASKLOOP_DEFINED -D NUM_TASKS=$(NUMTASKS)
ONLINECORESFLAG=-D ONLINECORES=$(ONLINE_CORES)


# (PARFOR_STATIC PARFOR_GUIDED PARFOR_DYNAMIC TASKLOOP)
SRC=src
UTIL=util
UTILS=$(UTIL)/*.c
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
	$(CC) -D $(CAPABILITY) $(ONLINECORESFLAG) $(INC) $(LDFLAGS) $(TASKFL) $(UTILS) $@ $(EM) -o $(BIN)/$(basename $(notdir $@))_task
else
$(SRCS): 
endif
	$(CC) -D $(CAPABILITY) $(ONLINECORESFLAG) $(INC)  $(LDFLAGS)   $(UTILS) $@  $(EM) -o $(BIN)/$(basename $(notdir $@))_static  
	$(CC) -D $(CAPABILITY) $(ONLINECORESFLAG) $(INC) $(LDFLAGS) $(DYNAMICFL)  $(UTILS) $@ $(EM) -o $(BIN)/$(basename $(notdir $@))_dynamic 
	$(CC) -D $(CAPABILITY) $(ONLINECORESFLAG) $(INC) $(LDFLAGS) $(GUIDEDFL)  $(UTILS) $@ $(EM) -o $(BIN)/$(basename $(notdir $@))_guided


preprocess:
	$(CC) -D $(CAPABILITY) $(INCWITOUTBIGLITTLE)  $(LDFLAGS) -D PAR_CHUNKSIZE=1024  -o $(BIN)/preprocess  $(UTIL)/*.c $(SRC)/preprocess.c
	$(CPP) -D $(CAPABILITY) $(INCWITOUTBIGLITTLE)  $(LDFLAGS) -D PAR_CHUNKSIZE=1024  -o $(BIN)/graphgenerator  $(UTIL)/*.c  $(SRC)/graphProperty.cpp  $(SRC)/graphgenerator.cpp -lsprng  ${SPRNG_BIN[@]/#/'/usr/local/bin/'}


dynamicwithchunkSize:
	$(CC) -D $(CAPABILITY) $(ONLINECORESFLAG) $(INC) $(LDFLAGS) $(DYNAMICFL) $(TARGET) $(EM) -o $(BIN)/$(basename $(notdir $(TARGET)))_dynamic_$(CHUNKSIZE) 

dynamicwithchunkSizedebug:
	$(CC) -g -D $(CAPABILITY) $(ONLINECORESFLAG) $(INC) $(LDFLAGS) $(DYNAMICFL) $(TARGET) $(EM) -o $(BIN)/$(basename $(notdir $(TARGET)))_dynamic_$(CHUNKSIZE) 




fromprecompiled: precompile
	$(CC) $(LDFLAG) -o $(BIN)/$(basename $(notdir $<)) $<

ifeq ( $(TASKLOOP_DEFINED), yes )
precompile: src/%.c
	$(CC) $(PRECOMPILE) $(INC) -D $(CAPABILITY) $(ONLINECORESFLAG) $(TASKFL) -o  $(PRE)/$(basename $(notdir $<))_task.i  $(UTILS) $<
else
precompile: src/%.c
endif
	$(CC) $(PRECOMPILE) $(INC) -D $(CAPABILITY) $(ONLINECORESFLAG) -o  $(PRE)/$(basename $(notdir $<))_static.i  $(UTILS) $<
	$(CC) $(PRECOMPILE) $(INC) -D $(CAPABILITY) $(ONLINECORESFLAG) $(DYNAMICFL) -o  $(PRE)/$(basename $(notdir $<))_dynamic.i $(UTILS) $<
	$(CC) $(PRECOMPILE) $(INC) -D $(CAPABILITY) $(ONLINECORESFLAG) $(GUIDEDFL) -o  $(PRE)/$(basename $(notdir $<))_guided.i  $(UTILS) $<

clean:
	rm -f $(OBJ)/*
	rm -f $(BIN)/*
	rm -f $(PRE)/*
	rm -f $(DEBUG)/*


