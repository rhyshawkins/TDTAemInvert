

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

syntheticstudy.obs : syntheticstudy.image standard16x16.path $(LMNOISE) $(HMNOISE)
	../mksyntheticobservations -i syntheticstudy.image -I standard16x16.path \
	-S ../stm/SkytemLM-BHMAR.stm -S ../stm/SkytemHM-BHMAR.stm \
	-N $(LMNOISE) -N $(HMNOISE) \
	-o $@ -O syntheticstudy.true

syntheticstudy.model : syntheticstudy.image
	../analysemodel -i syntheticstudy.rawimage -d 4 -l 4 -w 4 -W 4 -i syntheticstudy.rawimage -t 0.025 -T $@ -L

syntheticstudy.like : syntheticstudy.obs syntheticstudy.model
	../modellikelihood -i syntheticstudy.model -o syntheticstudy.obs \
	-s ../stm/SkytemLM-BHMAR.stm -s ../stm/SkytemHM-BHMAR.stm \
	-D 200.0 -d 4 -l 4 \
	-w 4 -W 4 \
	-H $(LMNOISE) -H $(HMNOISE) > $@

synthetic_smooth_threshold_0.02 : syntheticstudy.rawimage
	../analysemodel -i $< -t 0.02 -T $@ -d 4 -l 4 -w 4 -W 4




clean :
	rm -f $(TARGETS)
