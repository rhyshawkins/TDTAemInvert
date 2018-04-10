

LMNOISE = ../noise_models/brodienoiseLM.txt
HMNOISE = ../noise_models/brodienoiseHM.txt

TARGETS = synthetic_smooth_threshold_0.02 \
	standard16x16.path \
	syntheticstudy.image \
	syntheticstudy.model \
	syntheticstudy.obs

all : $(TARGETS)

#
# Smooth image from simple layered model
#
synthetic_smooth : synthetic_pixels
	python2 smoothimage.py -i $< -o $@ -s 1.5

#
# Convert image to required format for synthetic observation generation
#
syntheticstudy.image : synthetic_smooth
	python2 convertimage.py -i synthetic_smooth -o $@ -r syntheticstudy.rawimage -d 200.0 --min 0.05 --max 0.20

#
# Generate synthetic flight path
#
standard16x16.path :
	../mksyntheticflightpath -o $@ \
	-N 16 \
	-e 30.0 -E 2.5 \
	-p -3.2 -P 2.0 \
	-r 0.0 -R 2.0 \
	-x -12.5 -X 0.1 \
	-z 2.0 -Z 0.1 \

#
# Generate synthetic observations
#
syntheticstudy.obs : syntheticstudy.image standard16x16.path $(LMNOISE) $(HMNOISE)
	../mksyntheticobservations -i syntheticstudy.image -I standard16x16.path \
	-S ../stm/SkytemLM-BHMAR.stm -S ../stm/SkytemHM-BHMAR.stm \
	-N $(LMNOISE) -N $(HMNOISE) \
	-o $@ -O syntheticstudy.true

#
# Construct approximate model (wavelet thresholded)
#
syntheticstudy.model : syntheticstudy.image
	../analysemodel -i syntheticstudy.rawimage -d 4 -l 4 -w 4 -W 4 -i syntheticstudy.rawimage -t 0.10 -T $@ -L


#
# Check likelihood of approximate model (likelihood should be around 340 ((19 + 21) * 16/2))
#
syntheticstudy.like : syntheticstudy.obs syntheticstudy.model
	../modellikelihood -i syntheticstudy.model -o syntheticstudy.obs \
	-s ../stm/SkytemLM-BHMAR.stm -s ../stm/SkytemHM-BHMAR.stm \
	-D 200.0 -d 4 -l 4 \
	-w 4 -W 4 \
	-H $(LMNOISE) -H $(HMNOISE) > $@


results/ch.dat : syntheticstudy.obs
	mkdir -p results
	../aeminvert -i syntheticstudy.obs \
	-I syntheticstudy.model \
	-s ../stm/SkytemLM-BHMAR.stm -s ../stm/SkytemHM-BHMAR.stm \
	-M prior_laplacian.txt \
	-d 4 -l 4 -w 4 -W 4 \
	-H $(LMNOISE) -H $(HMNOISE) \
	-o results/ \
	-t 1000 -v 10 \
	-D 200.0

results/mean.txt : results/ch.dat
	../postprocess_mean -d 4 -l 4 \
	-i results/ch.dat \
	-o results/mean.txt \
	-D results/stddev.txt \
	-w 4 -W 4 \
	-t 5 -s 500

results_pt/ch.dat-000 : syntheticstudy.obs
	mkdir -p results_pt
	mpirun -np 8 ../aeminvert_pt -i syntheticstudy.obs \
	-s ../stm/SkytemLM-BHMAR.stm -s ../stm/SkytemHM-BHMAR.stm \
	-M prior_laplacian.txt \
	-d 4 -l 4 -w 4 -W 4 \
	-H $(LMNOISE) -H $(HMNOISE) \
	-L 0.1 -p 0.1 \
	-c 2 -T 4 -m 5.0 \
	-o results_pt/ \
	-t 10000 -v 100 \
	-D 200.0

results_pt_cont/ch.dat-000 : syntheticstudy.obs results_pt/ch.dat-000
	mkdir -p results_pt_cont
	mpirun -np 8 ../aeminvert_pt -i syntheticstudy.obs \
	-I results_pt/ \
	-s ../stm/SkytemLM-BHMAR.stm -s ../stm/SkytemHM-BHMAR.stm \
	-M prior_laplacian.txt \
	-d 4 -l 4 -w 4 -W 4 \
	-H $(LMNOISE) -H $(HMNOISE) \
	-L 0.1 -p 0.1 \
	-c 2 -T 4 -m 5.0 \
	-o results_pt_cont/ \
	-t 10000 -v 100 \
	-D 200.0

results_pt/mean.txt : results_pt/ch.dat-000
	mpirun -np 2 ../postprocess_mean_mpi -d 4 -l 4 \
	-i results_pt/ch.dat \
	-o results_pt/mean.txt \
	-D results_pt/stddev.txt \
	-w 4 -W 4 \
	-t 5 -s 5000

results_pt_cont/mean.txt : results_pt_cont/ch.dat-000
	mpirun -np 2 ../postprocess_mean_mpi -d 4 -l 4 \
	-i results_pt_cont/ch.dat \
	-o results_pt_cont/mean.txt \
	-D results_pt_cont/stddev.txt \
	-w 4 -W 4 \
	-t 5 -s 5000


results_pt_cont/like.txt : results_pt_cont/ch.dat-000
	../postprocess_likelihood \
	-i results_pt_cont/ch.dat-000 \
	-o results_pt_cont/like.txt \
	-H results_pt_cont/hierarchical.txt \
	-P results_pt_cont/prior.txt \
	-t 5

clean :
	rm -f $(TARGETS)
