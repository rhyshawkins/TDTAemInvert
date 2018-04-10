
import argparse

import numpy
import scipy.ndimage

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', type = str, required = True, help = 'Input image')
    parser.add_argument('-o', '--output', type = str, required = True, help = 'Output image')

    parser.add_argument('-s', '--sigma', type = float, default = 1.0, help = 'Sigma of Gaussian filter')
    
    args = parser.parse_args()

    img = numpy.loadtxt(args.input)

    fimg = scipy.ndimage.gaussian_filter(img, args.sigma)

    numpy.savetxt(args.output, fimg)

    
    
