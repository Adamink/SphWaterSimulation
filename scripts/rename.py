import os
from glob import glob
fd = './output/c0.5p0.05k10.0s0.6/'
for pth in glob(fd + '*'):
    file_name = os.path.basename(pth)
    index = int(file_name.split('.')[0])
    if index < 1000:
        new_file_name = '{:04d}.obj'.format(index)
        new_pth = os.path.join(fd, new_file_name)
        os.system(f'mv {pth} {new_pth}')
