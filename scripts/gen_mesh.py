import os
from datetime import datetime


class MeshGenerator:
    def __init__(self, debug=False, parent_fd='./output/', sub_fd='liquid',
                 cube_size=0.5, particle_radius=0.05, kernel_radius=10.0, sufrace_threshold=0.6):

        setting_str = f"c{cube_size}p{particle_radius}k{kernel_radius}s{sufrace_threshold}_5"
        output_fd = os.path.join(parent_fd, "meshed_liquid", setting_str)
        self.log_file = os.path.join(parent_fd, setting_str + '.log')

        if not os.path.exists(output_fd):
            os.makedirs(output_fd)

        param_str = f" --cube-size={cube_size} --particle-radius={particle_radius}" + \
            f" --kernel-radius={kernel_radius} --surface-threshold={sufrace_threshold}"

        if debug:
            frameRange = range(10, 20)
        else:
            frameRange = range(1030, 5000, 5)

        for frame in frameRange:
            file_pth = os.path.join(
                parent_fd, sub_fd, "{:04d}.ply".format(frame))
            index = os.path.basename(file_pth).split('.')[0]
            vtk_pth = os.path.join(parent_fd, index + '.vtk')
            output_pth = os.path.join(output_fd, index + ".obj")

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
    MeshGenerator()
