# SphWaterSimulation
This project is a course work of Physically-Based Simulation in Computer Graphics, ETH Zurich, 2020 Fall.

The is a collaborative project from [Ge Cao](https://github.com/GeCao), [Xiao Wu](https://github.com/Adamink), [Mingyang Song](https://github.com/FrauSong).

The demo video is available at Youtube.


## Install & Run
The installation has only been tested on Ubuntu for now.

The project is built with [cmake](https://cmake.org/). Before installation, please make sure cmake is installed.

To build the project, run the following command in the project folder.
```shell
mkdir build
cd build
cmake ..
make
```
Then, type `./SphWaterSimulation` to run the simulation.

## Rendering
For further rendering, we use [splashsurf](https://github.com/w1th0utnam3/splashsurf) to reconstruct the liquid surface, and use [blender](https://www.blender.org/) to render the scene. When importing sequence of .obj files into blender, please refer to plugin [stop motion obj](https://github.com/neverhood311/Stop-motion-OBJ) for blender.

To install relevant packages, you can refer to `scripts/install_utilities.sh`

We also upload our blender file for future reference. [Link](https://drive.google.com/drive/folders/1yZUP7o5rQNcQyGSNhPJRxcLfqKRuhjfM?usp=sharing).


