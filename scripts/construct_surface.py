import os
from datetime import datetime


class SurfaceConstructor:
    def __init__(self, frameRange=range(1, 201),
                 parent_fd='./output/', in_fd='liquid_ply', out_fd='meshed_liquid',
                 cube_size=0.5, particle_radius=0.05, kernel_radius=10.0, surface_threshold=0.6):

        setting_str = f"c{cube_size}p{particle_radius}k{kernel_radius}s{surface_threshold}"
        sub_fd = os.path.join(parent_fd, out_fd, setting_str)
        self.log_file = os.path.join(parent_fd, out_fd, setting_str + '.log')

        if not os.path.exists(sub_fd):
            os.makedirs(sub_fd)

        param_str = f" --cube-size={cube_size} --particle-radius={particle_radius}" + \
            f" --kernel-radius={kernel_radius} --surface-threshold={surface_threshold}"

        for frame in frameRange:
            file_pth = os.path.join(
                parent_fd, in_fd, "{:04d}.ply".format(frame))
            index = os.path.basename(file_pth).split('.')[0]
            vtk_pth = os.path.join(sub_fd, index + '.vtk')
            output_pth = os.path.join(sub_fd, index + ".obj")

            self.log("Processing frame {:04d}...".format(frame))
            if not os.path.exists(file_pth):
                self.log(
                    "Stop with no input file {:s}.".format(file_pth))
                break
            if os.path.exists(output_pth):
                self.log("Ignore existed frame {:04d}.".format(frame))
                continue

            os.system("splashsurf " + file_pth + " -o " + vtk_pth + param_str)
            if not os.path.exists(vtk_pth):
                self.log(
                    "Surface construction fails on frame {:04d}".format(frame))
                os.system("rm -rf " + vtk_pth)
                continue

            os.system("meshio-convert " + vtk_pth + " " + output_pth)
            if not os.path.exists(output_pth):
                self.log(
                    "Converting from .vtk to .obj fails on frame {:04d}".format(frame))

            os.system("rm -rf " + vtk_pth)
        self.log("Surface construction finished.")

    def log(self, message: str):
        with open(self.log_file, 'a+') as f:
            timestamp = datetime.now().strftime("[%H:%M:%S] ")
            print(timestamp + message, file=f)


if __name__ == '__main__':
    SurfaceConstructor()
