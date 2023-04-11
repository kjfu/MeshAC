<!--
 * @Author: Kejie Fu
 * @Date: 2023-03-11 23:20:09
 * @LastEditTime: 2023-04-03 16:51:45
 * @LastEditors: Kejie Fu
 * @Description: 
 * @FilePath: /MeshAC/README.md
-->

# MeshAC

A 3D Mesh Generation and Adaptation Package for Multiscale (atomistic-to-continuum) Coupling Simulation for Materials Defects. Written by Kejie Fu, Mingjie Liao and Yangshuai Wang.

![two voids](./data/two_voids.jpg)

## Getting Started

You can install MeshAC with ...
```
...
```

## Overview

We now summarize the main components of the library. 

1. 

2. 

3. 

The following functionals are currently supported:
- [`...`](...) ....

### A/C coupling method in Julia

Atomistic-to-continuum coupling method in Julia. The multi-scale model is coupled by BGFC method.

The atomistic model is based on the Julia package JuLIP.

Please do remember to modify the path of mesher3d in AtC constructor (AtC.jl).

Example:

```julia
...
```

## Development

MeshAC is under active development. Please don't hesitate to open feature requests to help us guide development. We more than welcome contributions!

## Publications

MeshAC has been used in the following publications.

1. [Adaptive Multigrid Strategy for Geometry Optimization of Large-Scale Three Dimensional Molecular Mechanics (J. Comp. Phys. 2023)](https://www.sciencedirect.com/science/article/pii/S0021999123002085)<br> K. Fu, M. Liao, Y. Wang, J. Chen and L. Zhang

## Citation

If you use the codes in a publication, please cite the repo using the .bib,

```
@inproceedings{...,
 author = {...},
 booktitle = {...},
 publisher = {...},
 title = {...},
 url = {...},
 volume = {...},
 year = {...}
}
```
