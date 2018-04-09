
TDTBASE=../TDTbase

#
# TODO Remove GAAEM from PhD repo
#
GAAEMBASE=/home/rhys/PhD/PhDS/libraries

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

EXTRADIST = libraries/Makefile \
	libraries/ga-aem.tar.gz \
	libraries/hnk.tar.gz \
	libraries/log.tar.gz \
	libraries/oset.tar.gz \
	libraries/sphericalwavelet.tar.gz \
	libraries/tracking.tar.gz \
	libraries/wavelet.tar.gz \
	libraries/wavetree.tar.gz \
	synthetic/Makefile \
	synthetictests/Makefile \
	syntheticstudy/Makefile \
	syntheticstudy/synthetic_smooth \
	syntheticstudy/synthetic_pixels \
	syntheticstudy/convertimage.py \
	syntheticstudy/syntheticstudy.prior \
	syntheticstudy/syntheticstudy.prior.dw \
	data/23330.asc \
	data/23330.obs \
	data/23330_32_0.obs \
	data/23330_32_1.obs \
	data/23330_32_2.obs \
	data/23330_32_3.obs \
	data/23330_128.obs \
	data/23330_256.obs \
	data/23330_512.obs \
	data/23330_threshold_0.25.model \
	data/23330_128_threshold_0.25.model \
	data/23330_256_threshold_0.25.model \
	data/23330_512_threshold_0.25.model \
	data/prior.txt \
	data/prior_coarse.txt \
	data/prior_final.txt \
	data/prior_fine.txt \
	data/prior32_final.txt \
	data/prior32_laplacian.txt \
	data/hypnoiseLM.txt \
	data/hypnoiseHM.txt \
	data/brodienoiseLM.txt \
	data/brodienoiseHM.txt \
	data/cov32_0_1HM.txt \
	data/cov32_0_1LM.txt \
	data/cov32_0_1HM_1.txt \
	data/cov32_0_1LM_1.txt \
	data/cov32_2_1HM.txt \
	data/cov32_2_1LM.txt \
	datatests/Makefile \
	stm/SkytemHM-BHMAR.stm \
	stm/SkytemLM-BHMAR.stm \
	stm/Skytem-HM.stm \
	stm/Skytem-LM.stm \
	stm/Tempest-standard.stm \
	scripts/convert.py \
	terrawulf/pbs_23330_init.sh \
	terrawulf/pbs_23330_256_init.sh \
	terrawulf/pbs_23330_256_init_nocache.sh \
	terrawulf/pbs_23330_256_init_fromzero_nocache.sh \
	terrawulf/pbs_23330_mpi_init.sh \
	terrawulf/pbs_23330_8_4_6_init.sh \
	terrawulf/pbs_23330_8_4_12_init.sh \
	terrawulf/pbs_23330_8_4_12_cont.sh \
	terrawulf/pbs_23330_32_6_5_4_init.sh \
	terrawulf/pbs_23330_32_6_5_4_cont.sh \
	terrawulf/pbs_23330_32_6_5_4_init_cov1.sh \
	terrawulf/pbs_23330_32_6_5_4_cont_cov1.sh \
	terrawulf/pbs_23330_32_6_10_6_init.sh \
	terrawulf/pbs_23330_32_6_10_6_cont.sh \
	terrawulf/pbs_23330_4_2_24_init.sh \
	terrawulf/pbs_syntheticstudy_4_3_4_init.sh \
	terrawulf/pbs_syntheticstudy_4_3_4_cont.sh \
	terrawulf/pbs_syntheticstudy_4_3_4_init_nonh.sh \
	terrawulf/pbs_syntheticstudy_4_3_4_cont_nonh.sh \
	terrawulf/pbs_syntheticstudy_4_4_4_T2_init.sh \
	terrawulf/pbs_syntheticstudy_4_4_4_T2_cont.sh \
	terrawulf/pbs_syntheticstudy_4_4_4_T2dw_init.sh \
	terrawulf/pbs_syntheticstudy_4_4_4_T2dw_cont.sh \
	terrawulf/8_4_12.sh \
	terrawulf/4_4_4_T2.sh \
	terrawulf/4_4_4_T2dw.sh \
	terrawulf/pbs_constant8x8standard_init.sh \
	terrawulf/pbs_23330_brodie_8_4_12_init.sh \
	terrawulf/pbs_23330_brodie_8_4_12_cont.sh \
	terrawulf/pbs_23330_128_8_4_12_init.sh \
	terrawulf/pbs_23330_128_8_4_12_cont.sh \
	terrawulf/23330_32_6_10_6_0.sh \
	terrawulf/23330_32_6_5_4_0.sh \
	terrawulf/23330_32_6_5_4_1.sh \
	terrawulf/23330_32_6_5_4_2.sh \
	terrawulf/23330_32_6_5_4_3.sh \
	terrawulf/pbs_23330_32_6_10_6_init_haar.sh \
	terrawulf/pbs_23330_32_6_10_6_cont_haar.sh \
	raijin/pbs_23330_init.sh \
	raijin/pbs_23330_8_4_16_init.sh \
	raijin/pbs_23330_8_4_16_cont.sh \
	raijin/pbs_23330_32_8_4_4_init.sh \
	raijin/pbs_23330_32_8_4_4_cont.sh \
	raijin/pbs_23330_32_8_8_4_init.sh \
	raijin/pbs_23330_32_8_8_4_cont.sh \
	raijin/pbs_23330_32_8_8_4_hp_init.sh \
	raijin/pbs_23330_32_8_8_4_hp_cont.sh \
	raijin/pbs_23330_32_8_8_4_hp_cov_init.sh \
	raijin/pbs_23330_32_8_8_4_hp_cov_cont.sh \
	raijin/pbs_23330_32_8_8_4_cov_init.sh \
	raijin/pbs_23330_32_8_8_4_cov_cont.sh \
	raijin/pbs_23330_128_8_4_16_init.sh \
	raijin/pbs_23330_128_8_4_16_cont.sh \
	raijin/pbs_23330_128_8_4_16_init_hierarchical.sh \
	raijin/pbs_23330_128_8_4_16_cont_hierarchical.sh \
	raijin/pbs_23330_128_8_4_16_cont_hierarchical_sample.sh \
	raijin/pbs_23330_8_4_16_test.sh \
	raijin/pbs_23330_32_8_32_4_cov_cont.sh \
	raijin/pbs_23330_32_8_32_4_cov_init.sh \
	raijin/23330_32_0.sh \
	raijin/23330_32_1.sh \
	raijin/23330_32_2.sh \
	raijin/23330_32_3.sh \
	raijin/23330_32_cov_0.sh \
	raijin/23330_32_32_0.sh \
	raijin/cov32_0_HP_LM.txt \
	raijin/cov32_0_HP_HM.txt 


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
DIR = AemInvert
TGZ = $(DIR).tar.gz

dist :
	make -C libraries dist
	mkdir -p $(DIR)
	echo $(DATE) > $(DIR)/Version
	for f in Makefile $(SRCS) $(EXTRADIST); do \
	    $(INSTALL) $(INSTALLFLAGS) $$f $(DIR)/$$f ; \
	done
	tar -czf $(TGZ) $(DIR)/*
	rm -rf $(DIR)

clean :
	rm -f $(TARGETS) *.o




