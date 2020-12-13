import os
from glob import glob

input_fd = './output/meshed_liquid/c0.5p0.05k10.0s0.6_5/'
output_fd = './output/meshed_liquid/c0.5p0.05k10.0s0.6_5div5/'

if not os.path.exists(output_fd):
    os.makedirs(output_fd)

for pth in glob(os.path.join(input_fd, '*')):
    basename = os.path.basename(pth)
    index = int(basename.split('.')[0])
    postfix = basename.split('.')[-1]
    if index % 5 == 0 and index <= 1000:
        newindex = int(index / 5)
        newbasename = '{:04d}.'.format(newindex) + postfix
        # output_pth = os.path.join(output_fd, basename)
        output_pth = os.path.join(output_fd, newbasename)
        os.system(f'cp {pth} {output_pth}')