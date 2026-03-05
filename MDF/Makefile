#Super duper Makefile from Hans
CXX=g++
LINKER=g++

ROOTCONFIG := root-config

CFLAGS := $(shell $(ROOTCONFIG) --cflags)
CFLAGS += -I${FAIRROOTPATH}/include
CFLAGS += -I$(SIMPATH)/include
CFLAGS += -I$(ROOT_INCLUDE_PATH)
CFLAGS += -I$(ROOT_INCLUDE_DIR)
CFLAGS += -I$(VMCWORKDIR)
CFLAGS += -I$(VMCWORKDIR)/tracking
CFLAGS += -I$(VMCWORKDIR)/field
CFLAGS += -I$(VMCWORKDIR)/r3bsource/base
CFLAGS += -I$(VMCWORKDIR)/r3bbase
CFLAGS += -I$(UCESB_DIR)/hbook
CFLAGS += --std=c++17 -g -O0 -fexceptions
CFLAGS += -g -O0 -fexceptions
#CFLAGS += -Werror
#CFLAGS += -g -Wall -W -Wconversion -Wshadow -Wcast-qual -Wwrite-strings

LDFLAGS := $(shell $(ROOTCONFIG) --ldflags)
LDFLAGS := $(shell $(ROOTCONFIG) --libs)
LDFLAGS += -L$(ROOT_LIBRARY_DIR) -L$(FAIRROOTPATH)/lib
LDFLAGS += -g
LDFLAGS += -L$(VMCWORKDIR)/../build/lib -lR3BBase -lField -lR3BTracking
LDFLAGS += -L$(FAIRROOTPATH)/lib -lBase -lParBase -lFairTools

INCLUDEDIR=include
DIR_INC=-I$(INCLUDEDIR)
EXEC=MDF
OBJ=MDF.o run_MDF.o
HEADERS= definitions.hh	\
	 libs.hh

DEPS = $(patsubst %,$(INCLUDEDIR)/%,$(HEADERS))

all: $(EXEC)

$(EXEC): $(OBJ)
	${LINKER} -o $@ $^ ${LDFLAGS} ${CFLAGS} ${DIR_INC}
	echo " COMP $(EXEC)"

%.o : %.cpp $(DEPS)
	$(MAKEDEPEND)
	${CXX} ${CFLAGS} $(DIR_INC) -c $< -o $@
	echo "	CXX $@"


clean:
	rm -f *.o $(EXEC)
	rm -f MDF
