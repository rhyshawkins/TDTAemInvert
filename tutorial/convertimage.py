
import numpy
import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', type = str, required = True, help = 'Input raw image')
    parser.add_argument('-o', '--output', type = str, required = True, help = 'Output aem image')
    parser.add_argument('-r', '--raw-output', type = str, default = None, help = 'Output raw image')

    parser.add_argument('-d', '--depth', type = float, default = 200.0, help = 'Depth to halfspace')
    parser.add_argument('--min', type = float, default = 0.05, help = 'Min value')
    parser.add_argument('--max', type = float, default = 0.20, help = 'Max value')

    args = parser.parse_args()

    m = numpy.loadtxt(args.input)
    rows, cols = m.shape

    minv = numpy.min(m)
    maxv = numpy.max(m)

    #
    # Rescale Image
    #
    m = (m - minv)/(maxv - minv) * (args.max - args.min) + args.min

    print rows, cols
    
    f = open(args.output, 'w')

    f.write('%d %d %f\n' % (rows, cols, args.depth))

    for row in range(rows):

        f.write(' '.join(map(lambda x : '%15.9f ' % x, m[row,:])) + '\n')

    f.close()

    if args.raw_output:
        numpy.savetxt(args.raw_output, m)
        
