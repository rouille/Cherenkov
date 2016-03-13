############
# makefile # 
############

SYSTEM=$(shell uname)
CC = gcc 
CXX = g++
OPTFLAGS = -O2
DBGFLAGS = -g3
WFLAGS = -D__USE_FIXED_PROTOTYPES__ -Wall
OBJ = ./

.PHONY: clean
#------------------- defs ------------------------------

INCDIR = $(shell root-config --cflags)

LIBDIR = $(shell root-config --glibs)

#-------------------------------------------------------

#------- alias -----------------------------------------
libobjs = \
          common.o \
          atmosphere.o \
          cherenkov.o \
          conversion.o \
          shower.o 

execs = \
        example_atmosphere.exe \
        example_shower.exe \
        example_spectra.exe \
        example_cherenkov.exe 


exeobjs = $(patsubst %.exe,%.o,$(execs))

HEADERS = $(patsubst %.o,%.h,$(libobjs))

thelib = lib.a
#-------------------------------------------------------

#-------- rules ----------------------------------------
# rules for the library sources
$(OBJ)%.o:%.cc %.h
	$(COMPILE.cc) $(DBGFLAGS) $(OPTFLAGS) $(WFLAGS) $(INCDIR) -o $@ $<

# rules for the executable sources
$(OBJ)%.o:%.cc
	$(COMPILE.cc) $(DBGFLAGS) $(OPTFLAGS) $(WFLAGS) $(INCDIR) -o $@ $<
#-------------------------------------------------------

#------- targets ---------------------------------------
all : lib $(execs)
lib: $(thelib) 
clean :
	@echo "Deleting library objects, executables and associated objects."
	@/bin/rm -f $(thelib) $(execs) $(libobjs) $(exeobjs)
	@echo "Done"
#-------------------------------------------------------

#-------- specific rules -------------------------------
$(thelib): $(libobjs)
	@echo "Making "$(thelib)
	@ar r $@ $^
	@ranlib $@
	@echo "Done."
example_atmosphere.exe: example_atmosphere.o $(thelib)
	$(CXX) -o $@ $^ $(LIBDIR)
example_shower.exe: example_shower.o $(thelib)
	$(CXX) -o $@ $^ $(LIBDIR)
example_spectra.exe: example_spectra.o $(thelib)
	$(CXX) -o $@ $^ $(LIBDIR)
example_cherenkov.exe: example_cherenkov.o $(thelib)
	$(CXX) -o $@ $^ $(LIBDIR)
#-------------------------------------------------------

