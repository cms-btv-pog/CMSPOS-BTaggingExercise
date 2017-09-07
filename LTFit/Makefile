# Makefile

CC = g++ -g -Wno-write-strings

CFITINC    = -I../CFIT/
CFITLIB    = -L../CFIT/ -lCFIT

SOVERSION               = $(shell echo $(VERSION) | \
                               sed 's/\(.*\)\.\(.*\)\.\(.*\)/\1:\2:0/')
ROOTCONF                = root-config
ROOTCINT                = $(shell $(ROOTCONF) --prefix)/bin/rootcint
ROOTLIBS                = $(shell $(ROOTCONF) --libs)
ROOTGLIBS               = $(shell $(ROOTCONF) --glibs)
ROOTCFLAGS              = $(shell $(ROOTCONF) --cflags)
INCLUDES                = $(shell $(ROOTCONF) --incdir) $(ROOTCFLAGS)
CINTINCLUDES            = $(shell $(ROOTCONF) --incdir)
LIBS                    = $(ROOTLIBS)

LDLIBS  = $(ROOTLIBS) $(CFITLIB)

ADDINC = $(CFITINC)

SRC = test.cxx

INC = ../CFIT/libCFIT.so

OBJ = test.o

EXE = test

all:  	     		$(EXE)
#			@rm -f *.o

$(EXE):     	 	$(OBJ) $@

$(INC):
	cd ../CFIT/ && make && cd -

$(OBJ):     		$(SRC) $(INC)
		        $(CC) -c $(ROOTCFLAGS) $(ADDINC) $(SRC)

clean: 		
			@rm -f $(OBJ) $(EXE)
