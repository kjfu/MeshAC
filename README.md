<!--
 * @Author: Kejie Fu
 * @Date: 2023-03-11 23:20:09
 * @LastEditTime: 2023-04-23 23:22:17
 * @LastEditors: Kejie Fu
 * @Description: 
 * @FilePath: /MeshAC/README.md
-->

# MeshAC

A 3D Mesh Generation and Adaptation Package for Multiscale (atomistic-to-continuum) Coupling Simulation for Materials Defects. Written by Kejie Fu, Mingjie Liao and Yangshuai Wang.

![two voids](./data/two_voids.jpg)

## Getting Started
### Installation

You can install MeshAC with the following steps:
- Clone the repository into your local machine:

```bash
git clone https://github.com/kjfu/MeshAC.git
```

- Compile the code using cmake:
  
```bash
cd MeshAC
mkdir build
cd build
cmake ..
make
```
### Command Line Switches

MeshAC supports usage via command line. Here is an overview of all command line switches:

```
Usage: ./MeshAC [OPTIONS] input [output]

Postionals:
	input TEXT REQUIRED     Input a file of initial mesh(.mesh) for mesh generation or middle files
				(.mesh, .remesh and .value) for mesh adaptation.(string, required)
	output TEXT             Output a file of resulting mesh (.mesh) or a file of interpolation solutions
				(.value).

Options:
	-h			Print this help message and exit.
	-i TEXT REQUIRED        Input a file of initial mesh(.mesh) for mesh generation or middle files
				(.mesh, .remesh and .value) for mesh adaptation.(string, required)
	-o TEXT                 Output a file of resulting mesh (.mesh) or a file of interpolation solutions
				(.value).
	-s  FLOAT               Input the max sizing value for mesh generation.
	-hd                     Generate a mesh with edge dislocation.
	-r                      Refine an existing mesh adaptively.
	-rr                     Refine an existing mesh with edge dislocation adaptively.
```
### Examples
#### To generate 3d mesh from points（with 8 points as bounding points with label 1, and several atomistic points with label 0）
```
>> ./MeshAC -s 5 -i sample.mesh -o outmesh.mesh
```
#### To generate a 3d mesh from points with edge dislocation

```
>>./MeshAC -hd -i test3d.mesh -s 15 -o out3d.mesh
````

#### To remesh a 3d mesh adaptively
You must keep 3 files (*.mesh, *.remesh, *.value) in same path.
```
>>./MeshAC -r -i test3d -o out3d
```

#### To remesh a 3d mesh with edge dislocation adaptively
You must keep 3 files(*.mesh, *.remesh, and *.value) in same path 
```
>>./MeshAC -rr -i test3d -o out3d
```

### Label Meanings in .mesh file

#### Labels for nodes

| Label | Significance |Tip|
|:------|:-------|:-----|
|0|Nodes inside the atomic area||
|1|Nodes on the border of the continuous area||
|2|Nodes on the border of the atomic area||
|3|Nodes between the border of the continuous area and the border of the atomic area||
#### Labels for tets

| Label | Significance |Tip|
|:------|:-------|:-----|
|0|Tets of the atomic area||
|1|Tets of the continuous area||


<!-- ## Overview

We now summarize the main components of the library. 

1. Mesh generation

2. Mesh adaptation

3. ...

The following functionals are currently supported:
- [`...`](...) ....

From Mesher3DForSJTU package ... -->

## A/C coupling method in Julia

Atomistic-to-continuum coupling method in Julia. The current implementation is based on the [BGFC (Blended Ghost Force Correction)](https://epubs.siam.org/doi/10.1137/15M1020241) method.

The atomistic computations involved are heavily depends on the pure Julia package [JuLIP](https://github.com/JuliaMolSim/JuLIP.jl) (Julia Library for Interatomic Potentials).

1. Install Julia 1.10.1 from [here](https://julialang.org/downloads/); 
2. Install necessary registry and then open a new activate environment:
```
] registry add https://github.com/JuliaRegistries/General"; 
] registry add https://github.com/ACEsuit/ACEregistry"; 
] activate .
```
3. Run the following Julia command to install the required packages:
```
] add JuLIP, DelimitedFiles, Printf, NeighbourLists, QHull, Optim, LineSearches, SparseArrays, Isaac, PyCall
```

Tips:
1. Use ```Pkg.activate(".")``` to use a local project and set environment variable ```JULIA_PROJECT``` accordingly. 
2. Please do remember to modify the path of mesher3d in AtC constructor (AtC.jl).

### FIO

The module named `ACFIO` is used to cope with the geometrical operations.

The functionals of FIO are to read/write:   

- .mesh
- .remesh   
- .value   

and write .dump files for visualization.

Example:

```julia
fn = "a.mesh"
# X: atom positions and outer continuum points; Xtype: interface information
ACFIO.write_mesh(fn, X, Xtype)
# call `mesher3d` to build coupled mesh
ofn = "ac.mesh"
run(`$meshpath -s $h -i $fn -o $ofn`)
# X: nodes; T: mesh topology
X, T = ACFIO.read_mesh(ofn)
```

### AtC

Major struct contains geometrical and computational information.

To construct an AtC objective invokes:   
```
function AtC(Ra::Int64, bw::Int64, Lmsh, h; Rbuf=2, sp=:W, r0=rnn(:W), defects=:SingVac, meshpath="YOUR PATH")
```
Note that `meshpath` is the path of mesher toolkit that may differ from each devices.

Example:

```julia
# construct atomistic region with defects. R: radius; 
atdef = get_atdef(R)
# construct a/c coupling structure. h: mesh size; L: size of computational domain
atc0 = AtC(atdef, h, L)
```

## Development

MeshAC is under active development. Please don't hesitate to open feature requests to help us guide development. We more than welcome contributions!

## Publications

MeshAC has been used in the following publications.

1. [Adaptive Multigrid Strategy for Geometry Optimization of Large-Scale Three Dimensional Molecular Mechanics (J. Comp. Phys. 2023)](https://www.sciencedirect.com/science/article/pii/S0021999123002085)<br> K. Fu, M. Liao, Y. Wang, J. Chen and L. Zhang

## Citation

If you use the codes in a publication, please cite the repo using the .bib,

```
@inproceedings{fu2024meshac,
  title={MeshAC: A 3D Mesh Generation and Adaptation Package for Multiscale Coupling Methods},
  author={Fu, Kejie and Liao, Mingjie and Wang, Yangshuai and Chen, Jianjun and Zhang, Lei},
  journal={arXiv preprint arXiv:2402.09446},
  year={2024}
}
```

## Acknowledgements

We used several useful libraries in our implement and testing listed as follows. We would like to especially thank their authors for their great work and publishing the code.

- [Tetgen](http://www.tetgen.org)
- [Triangle](http://www.cs.cmu.edu/~quake/triangle.html)
- [AABB.cc](https://github.com/lohedges/aabbcc.git)
- [kdtree](https://github.com/jtsiomb/kdtree.git)
