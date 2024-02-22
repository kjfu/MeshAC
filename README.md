<!--
 * @Author: Kejie Fu
 * @Date: 2023-03-11 23:20:09
 * @LastEditTime: 2024-02-22 23:22:17
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
cp ./build/MeshAC ./MeshAC
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
] add JuLIP, DelimitedFiles, Printf, NeighbourLists, QHull, Optim, LineSearches, SparseArrays, Isaac, PyCall, Plots, ASE, ScatteredInterpolation
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

### AtC

Major struct contains geometry and computational information.

To construct an AtC objective invokes:   
```
function AtC(Ra::Int64, bw::Int64, Lmsh, h; Rbuf=2, sp=:W, r0=rnn(:W), defects=:SingVac, meshpath="YOUR PATH")
```
Note that `meshpath` is the path of mesher toolkit that may differ from each devices.

## Example
Step 1: Construct the atomistic model in atomistic region:
```julia
] activate .
include("./JuliaAC/dihole.jl")
# construct atomistic region with di-hole defects 
# L: repeated size (keep it sufficiently large); Lb: the inner ball size; Lw: the width of torus.
L = 20; Lb = 13.0; Lw = 16.0;
atdef = get_atdef(L, Lb, Lw);
```
Step 2: Write atomistic mesh:
```julia
Xat = positions(atdef)
cb = [-1 1 1 -1 -1 1 1 -1; -1 -1 1 1 -1 -1 1 1; -1 -1 -1 -1 1 1 1 1]
Xcb = (L + 20) .* cb
X = hcat(mat(Xat), Xcb)
Xtype = zeros(Int64, length(atdef))
append!(Xtype, ones(Int64, 8))
fn = "./data/dihole_A.mesh"
# X: atom positions and outer continuum points; Xtype: interface information
ACFIO.write_mesh(fn, X, Xtype)
```
Step 3: Build AtC mesh and read the coupled mesh in Julia:
```julia
# call `mesher3d` to build coupled mesh
ofn = "./data/dihole_AC.mesh"
# the density of nodes in continuum region
h = 10
run(`./MeshAC -s $h -i $fn -o $ofn`)
```
Step 4: Construct a/c coupling structure:
```julia
# X: nodes; T: mesh topology
X, T = ACFIO.read_mesh(ofn)
iBdry = findall(x->x==1.0, X[4,:])
XType = X[4,:]
nat = [findall(x->x==0.0, XType); findall(x->x==2.0, XType);]
Xat = deepcopy(X[1:3, nat])
set_positions!(atdef, Xat)
# alternative data
data = Dict{String, Real}()
U = zeros(3, size(X, 2))
∇U, volT, J = gradient(T, X, U)
wat = bulk(:W, cubic = true)
V0 = det(cell(wat))/2
# final a/c stucture
atc = AtC{eltype(X)}(atdef, V0, X[1:3, :], X[4, :], U, ∇U, T[1:4, :], T[5,:], 3, 2, J, wat, iBdry, data);
```
Step 5: Alternatively, one can call the following command directly:
```julia
] activate .
include("./JuliaAC/dihole.jl")
# collect model parameters
L = 20; Lb = 13.0; Lw = 16.0;
Lc = 20; h = 10;
input_name = "./data/double-voids.mesh"
output_name = "./data/double-voids-coupled.mesh"
# construct atomistic region with defects
atdef = get_atdef(L, Lb, Lw)
# construct a/c coupling structure
atc = AtC_di(atdef, h, L, Lc; fn = input_name, ofn = output_name);
```

## Development

MeshAC is under active development. Please don't hesitate to open feature requests to help us guide development. We more than welcome contributions!

## Publications

MeshAC has been used in the following publications.

1. [Adaptive Multigrid Strategy for Geometry Optimization of Large-Scale Three Dimensional Molecular Mechanics (J. Comp. Phys. 2023)](https://www.sciencedirect.com/science/article/pii/S0021999123002085)<br> K. Fu, M. Liao, Y. Wang, J. Chen and L. Zhang

## Citation

If you use the codes in a publication, please cite the repo using the .bib,

```
@misc{fu2024meshac,
      title={MeshAC: A 3D Mesh Generation and Adaptation Package for Multiscale Coupling Methods},
      author={Kejie Fu and Mingjie Liao and Yangshuai Wang and Jianjun Chen and Lei Zhang},
	  eprint={2402.09446},
      archivePrefix={arXiv},
      year={2024}
}
```

## Acknowledgements

We used several useful libraries in our implement and testing listed as follows. We would like to especially thank their authors for their great work and publishing the code.

- [Tetgen](http://www.tetgen.org)
- [Triangle](http://www.cs.cmu.edu/~quake/triangle.html)
- [AABB.cc](https://github.com/lohedges/aabbcc.git)
- [kdtree](https://github.com/jtsiomb/kdtree.git)
