# Adaptive module

## Adaptive Process

Solve $\to$ Estimate $\to$ Mark $\to$ Refine

### Solve

`BGFC/AtC` + `BGFC/Solve`

### Estimate

$\|\nabla u\|_{L^2(T)}$ for all $T\in\mathcal{T}$

### Mark

- Tetrahedral elements' indices.   
- Appending points' positions.

### Refine

Call `mesher3d -r`

## Call `mesher3d -r`

This refinement functional requires \*.mesh, \*.remash, and \*.value files in the same path.

- \*.mesh: original mesh files to be refined.
- \*.remesh: consists of two fields labeled with
	- Append\_points: vectors represent atomistic points to be appened adjacent to the interface.   
	- Refine_elements: indices of elements to be refined.
- \*.value: consists of two fileds labeled with (could be used individually)
	- scalar density
	- vector displacement