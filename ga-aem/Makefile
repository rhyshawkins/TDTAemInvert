
CXX = g++
CXXFLAGS = -c -g -Wall -std=c++11

CXXFLAGS += -O3

TARGETS = libga-aem.a

all : $(TARGETS)

OBJS = blocklanguage.o \
	file_utils.o \
	general_utils.o \
	geometry3d.o \
	le.o \
	tdemsystem.o

libga-aem.a : $(OBJS)
	$(AR) $(ARFLAGS) $@ $(OBJS)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -o $*.o $*.cpp

INSTALL = install
INSTALLFLAGS = -D

SRCS = blocklanguage.cpp \
	file_utils.cpp \
	general_utils.cpp \
	geometry3d.cpp \
	le.cpp \
	tdemsystem.cpp \
	blocklanguage.h \
	file_utils.h \
	general_constants.h \
	general_types.h \
	general_utils.h \
	geometry3d.h \
	le.h \
	tdemsystem.h \
	vector_utils.h

EXTRADIST = tests/Makefile \
	tests/Tempest-standard.stm \
	tests/test_single.cpp

DATE = $(shell date +"%Y%m%d%H%M")
DIR = ga-aem
TGZ = $(DIR).tar.gz

dist :
	mkdir -p $(DIR)
	echo $(DATE) > $(DIR)/Version
	for f in Makefile $(SRCS) $(EXTRADIST); do \
	    $(INSTALL) $(INSTALLFLAGS) $$f $(DIR)/$$f ; \
	done
	tar -czf $(TGZ) $(DIR)/*
	rm -rf $(DIR)


clean :
	rm -f $(TARGETS) *.o *~



