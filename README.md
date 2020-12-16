# SphWaterSimulation
This project is a course work of Physically-Based Simulation in Computer Graphics, ETH Zurich, 2020 Fall.

The is a collaborative project from Ge Cao, Xiao Wu, Mingyang Song.


## Install & Run
The installation has only been tested on Ubuntu.
To build the project, run the following command in the project folder.
```shell
mkdir build
cd build
cmake ..
make
```
Then, type `./SphWaterSimulation' to run the simulation.

For further rendering, we use [splashsurf](https://github.com/w1th0utnam3/splashsurf) to reconstruct liquid surface, and use [blender](https://www.blender.org/) to render the scene. While importing sequence of .obj files into blender, please refer to plugin [stop motion obj](https://github.com/neverhood311/Stop-motion-OBJ) for blender.



