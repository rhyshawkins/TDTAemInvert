
TDTBASE=../TDTbase

GAAEMBASE=.

INCLUDES = \
	$(shell mpicc -showme:compile) \
	-I$(TDTBASE)/log \
	-I$(TDTBASE)/hnk \
	-I$(TDTBASE)/tracking \
	-I$(GAAEMBASE)/ga-aem \
	-I$(TDTBASE)/wavetree \
	-I$(TDTBASE)/oset \
	-I$(TDTBASE)/wavelet 

EXTRA_LIBS = \
	-L$(TDTBASE)/hnk -lhnk \
	-L$(TDTBASE)/tracking -ltracking \
	-L$(TDTBASE)/log -llog \
	-L$(GAAEMBASE)/ga-aem -lga-aem \
	-lfftw3 \
	-L$(TDTBASE)/wavetree -lwavetree \
	-L$(TDTBASE)/oset -loset \
	-L$(TDTBASE)/wavelet -lwavelet

CXX = g++
CXXFLAGS = -c -g -Wall --std=c++11 -DOMPI_SKIP_MPICXX $(INCLUDES)

#ifeq ($(shell hostname),terrawulf)
CXXFLAGS += -O3
#endif

INSTALL = install
INSTALLFLAGS = -D

LIBS = $(EXTRA_LIBS) -lm $(shell gsl-config --libs) -lrt -lgmp 
MPI_LIBS = $(shell mpicxx -showme:link)

OBJS = aemexception.o \
	aemobservations.o \
	aemimage.o \
	aemutil.o \
	constants.o \
	hierarchical.o \
	hierarchicalmodel.o \
	hierarchicalprior.o \
	rng.o \
	logspace.o \
	global.o \
	global_pixel.o \
	birth.o \
	death.o \
	value.o \
	value_pixel.o \
	ptexchange.o \
	resample.o \
	chainhistory_pixel.o

MPI_OBJS = 

SRCS = Makefile \
	aemexception.cpp \
	aemimage.cpp \
	aeminvert.cpp \
	aeminvert_mpi.cpp \
	aeminvert_pixel.cpp \
	aeminvert_pt.cpp \
	aemobservations.cpp \
	aemutil.cpp \
	analysemodel.cpp \
	birth.cpp \
	chainhistory_pixel.cpp \
	computeresiduals.cpp \
	constants.cpp \
	death.cpp \
	global.cpp \
	global_pixel.cpp \
	hash.cpp \
	hierarchical.cpp \
	hierarchicalmodel.cpp \
	hierarchicalprior.cpp \
	logspace.cpp \
	mksyntheticflightpath.cpp \
	mksyntheticimage.cpp \
	mksyntheticobservations.cpp \
	modellikelihood.cpp \
	postprocess_coeff_marginal.cpp \
	postprocess_likelihood.cpp \
	postprocess_mean.cpp \
	postprocess_mean_mpi.cpp \
	postprocess_pixel_mean.cpp \
	postprocess_validate_likelihood.cpp \
	postprocess_acceptance.cpp \
	postprocess_coeff_history.cpp \
	postprocess_khistory.cpp \
	ptexchange.cpp \
	resample.cpp \
	rng.cpp \
	value.cpp \
	value_pixel.cpp \
	aemexception.hpp \
	aemimage.hpp \
	aemobservations.hpp \
	aemutil.hpp \
	birth.hpp \
	chainhistory_pixel.hpp \
	constants.hpp \
	death.hpp \
	global.hpp \
	global_pixel.hpp \
	hash.hpp \
	hierarchical.hpp \
	hierarchicalmodel.hpp \
	hierarchicalprior.hpp \
	logspace.hpp \
	ptexchange.hpp \
	resample.hpp \
	rng.hpp \
	value.hpp \
	value_pixel.hpp 

EXTRADIST = noise_models/brodienoiseHM.txt  \
	noise_models/brodienoiseLM.txt \
	stm/Skytem-HM.stm \
	stm/Skytem-LM.stm \
	stm/SkytemHM-BHMAR.stm \
	stm/SkytemLM-BHMAR.stm \
	stm/Tempest-standard.stm \
	doc/manual.tex \
	doc/bibliography.bib \
	doc/manual.pdf \
	tutorial/Makefile \
	tutorial/convertimage.py \
	tutorial/prior_laplacian.txt \
	tutorial/smoothimage.py \
	tutorial/synthetic_pixels

TARGETS = mksyntheticflightpath \
	mksyntheticimage \
	mksyntheticobservations \
	aeminvert \
	aeminvert_mpi \
	aeminvert_pt \
	aeminvert_pixel \
	postprocess_coeff_marginal \
	postprocess_likelihood \
	postprocess_mean \
	postprocess_mean_mpi \
	postprocess_pixel_mean \
	postprocess_validate_likelihood \
	postprocess_acceptance \
	postprocess_coeff_history \
	postprocess_khistory \
	analysemodel \
	modellikelihood \
	computeresiduals

all : $(TARGETS)

mksyntheticimage : mksyntheticimage.o $(OBJS)
	$(CXX) -o mksyntheticimage mksyntheticimage.o $(OBJS) $(LIBS) $(MPI_LIBS)

mksyntheticflightpath : mksyntheticflightpath.o $(OBJS)
	$(CXX) -o mksyntheticflightpath mksyntheticflightpath.o $(OBJS) $(LIBS) $(MPI_LIBS)

mksyntheticobservations : mksyntheticobservations.o $(OBJS)
	$(CXX) -o mksyntheticobservations mksyntheticobservations.o $(OBJS) $(LIBS) $(MPI_LIBS)

aeminvert : aeminvert.o $(OBJS)
	$(CXX) -o aeminvert aeminvert.o $(OBJS) $(LIBS) $(MPI_LIBS)

aeminvert_mpi : aeminvert_mpi.o $(OBJS)
	$(CXX) -o aeminvert_mpi aeminvert_mpi.o $(OBJS) $(LIBS) $(MPI_LIBS)

aeminvert_pt : aeminvert_pt.o $(OBJS)
	$(CXX) -o aeminvert_pt aeminvert_pt.o $(OBJS) $(LIBS) $(MPI_LIBS)

aeminvert_pixel : aeminvert_pixel.o $(OBJS)
	$(CXX) -o aeminvert_pixel aeminvert_pixel.o $(OBJS) $(LIBS) $(MPI_LIBS)

postprocess_coeff_marginal : postprocess_coeff_marginal.o $(OBJS)
	$(CXX) -o postprocess_coeff_marginal postprocess_coeff_marginal.o $(OBJS) $(LIBS) $(MPI_LIBS)

postprocess_likelihood : postprocess_likelihood.o $(OBJS)
	$(CXX) -o postprocess_likelihood postprocess_likelihood.o $(OBJS) $(LIBS) $(MPI_LIBS)

postprocess_mean : postprocess_mean.o $(OBJS)
	$(CXX) -o postprocess_mean postprocess_mean.o $(OBJS) $(LIBS) $(MPI_LIBS)

postprocess_mean_mpi : postprocess_mean_mpi.o $(OBJS)
	$(CXX) -o postprocess_mean_mpi postprocess_mean_mpi.o $(OBJS) $(LIBS) $(MPI_LIBS)

postprocess_pixel_mean : postprocess_pixel_mean.o $(OBJS)
	$(CXX) -o postprocess_pixel_mean postprocess_pixel_mean.o $(OBJS) $(LIBS) $(MPI_LIBS)

postprocess_validate_likelihood : postprocess_validate_likelihood.o $(OBJS)
	$(CXX) -o postprocess_validate_likelihood postprocess_validate_likelihood.o $(OBJS) $(LIBS) $(MPI_LIBS)

postprocess_acceptance : postprocess_acceptance.o $(OBJS)
	$(CXX) -o postprocess_acceptance postprocess_acceptance.o $(OBJS) $(LIBS) $(MPI_LIBS)

postprocess_coeff_history : postprocess_coeff_history.o $(OBJS)
	$(CXX) -o postprocess_coeff_history postprocess_coeff_history.o $(OBJS) $(LIBS) $(MPI_LIBS)

postprocess_khistory : postprocess_khistory.o $(OBJS)
	$(CXX) -o postprocess_khistory postprocess_khistory.o $(OBJS) $(LIBS) $(MPI_LIBS)

analysemodel : analysemodel.o $(OBJS)
	$(CXX) -o analysemodel analysemodel.o $(OBJS) $(LIBS) $(MPI_LIBS)

modellikelihood : modellikelihood.o $(OBJS)
	$(CXX) -o modellikelihood modellikelihood.o $(OBJS) $(LIBS) $(MPI_LIBS)

computeresiduals : computeresiduals.o $(OBJS)
	$(CXX) -o computeresiduals computeresiduals.o $(OBJS) $(LIBS) $(MPI_LIBS)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -o $*.o $*.cpp

DATE = $(shell date +"%Y%m%d%H%M")
DIR = TDTAemInvert
TGZ = $(DIR).tar.gz

doc/manual.pdf : doc/manual.tex
	cd doc && pdflatex manual.tex && bibtex manual && pdflatex manual.tex && pdflatex manual.tex

dist :
	mkdir -p $(DIR)
	echo $(DATE) > $(DIR)/Version
	for f in Makefile $(SRCS) $(EXTRADIST); do \
	    $(INSTALL) $(INSTALLFLAGS) $$f $(DIR)/$$f ; \
	done
	tar -czf $(TGZ) $(DIR)/*
	rm -rf $(DIR)

clean :
	rm -f $(TARGETS) *.o




