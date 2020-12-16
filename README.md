# SphWaterSimulation
![](./assets/demo.gif)

This repository is a course project of Physically-Based Simulation in Computer Graphics, ETH Zurich, Fall 2020.

The project is a collaborative work from [Ge Cao](https://github.com/GeCao), [Xiao Wu](https://github.com/Adamink), [Mingyang Song](https://github.com/FrauSong).

The demo video is available at [Youtube](https://youtu.be/wf_0VI8c-Rg).

Features:
- Particle-based fluid simulation with SPH method
- Iterative SESPH
- Coupling liquid simulation with rigid body (the wheel)
- Dambreak scene & Waterwheel scene
- Particles import & export in `ply` and `cfg` format
- Multiple ways of visualization using [OVITO](http://www.ovito.org/), [OpenGL/GLUT](https://www.opengl.org/resources/libraries/glut/) and [blender](https://www.blender.org/)
- Surface reconstruction using [splashsurf](https://github.com/w1th0utnam3/splashsurf)
- Multithreading for accelerating simulation using [OpenMP](https://www.openmp.org/)
- Rendered with GPU accelerated ray-tracing using [Cycles](https://www.cycles-renderer.org/) 

## Install & Run
The installation has only been tested on Ubuntu for now. The project is built with [CMake](https://cmake.org/). Before installation, please make sure CMake is installed. The project also depends on OpenGL/GLUT, please use the following command to install dependencies.
```shell
sudo apt-get install libglu1-mesa-dev freeglut3-dev mesa-common-dev libx11-dev libxi-dev
```

To build the project, run the following command in the project folder.
```shell
mkdir build
cd build
cmake ..
make
```
Then, type `./SphWaterSimulation` to run the simulation.


## Render
For further rendering, we use [splashsurf](https://github.com/w1th0utnam3/splashsurf) to reconstruct the liquid surface, and use [blender](https://www.blender.org/) to render the scene. When importing sequence of .obj files into blender, please refer to plugin [stop motion obj](https://github.com/neverhood311/Stop-motion-OBJ) for blender.

To install relevant packages, please refer to `scripts/install_utilities.sh`

After installation, use `python3 ./scripts/construct_surface.py` to construct liquid surface.

We also upload our blender file to Google Drive for future reference [Link](https://drive.google.com/drive/folders/1yZUP7o5rQNcQyGSNhPJRxcLfqKRuhjfM?usp=sharing).

## Config
The config file is `SphWaterSimulation/constants.h`. To disable visualization, change variable `IF_VISUALIZE` to `false`. To disable rigid body(the wheel), change `kUseRigidBody` to `false`.
