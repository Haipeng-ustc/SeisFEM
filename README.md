# **S**eisFEM

By Haipeng Li @ USTC

Contact: haipengl@mail.ustc.edu.cn

![](./Doc/logo/logo.png)

```bash
The SeisFEM is a 2D Finite Element code to simulate the elastic wave propagations. 
All the algorithm design is finished during my study at the University of Tulsa 
supervised by Prof. Jingyi Chen. The FEM code was rewritten into C during my later study.
If you have any problem, please contact me by email: haipengl@mail.ustc.edu.cn 
```



# Descriptions of each folder
## assemble 

Assemble the mass and stiffness matrices. All the matrices are stored in the csr format. (csr: compressed sparse row)
You can use different element types and here is the element type list:

	     I  ELEMENT_TYPE   Definition
		   -  ------------   ----------
	     1  T3              3 node linear triangle;
	     2  T6              6 node quadratic triangle;
	     3  T10            10 node cubic triangle.
	     4  Q4              4 node linear Lagrange/serendipity quadrilateral;
	     5  Q9              9 node quadratic Lagrange quadrilateral;
	     6  Q16            16 node cubic Lagrange quadrilateral;  
	     
	Note: T10 and Q16 element types require too many memories and take ages to finish the simulation.

## backup 

# Reference
Li, H., Chen, J., Zhao, Z., & Li, J. (2021). A multi-axial perfectly matched layer for finite-element
time-domain simulation of anisotropic elastic wave propagation. Journal of Seismic Exploration,
30, 173-200.

Old codes, just in case need to use them one day.

## Doc

Some tutorials on the Finite Element Method.

## Example

Provided examples. You can tell their simulation types by their names:

	  example1_first_order_triangular_mesh
	  example2_first_order_quadrilateral_mesh
	  example3_second_order_triangular_mesh
	  example4_second_order_quadrilateral_mesh
	  example5_use_exterior_mesh
	  example6_surface_wave
	  example7_topography_surface

## femlib: 

Finite Element Method library from: https://people.sc.fsu.edu/~jburkardt/f_src/f_src.html.

    Both Fortran and C libraries are included here. 
    They can be very helpful for the beginner of the Finite Element Method to write their own code.

## mesh: 

Perform the mesh of the computational domain using the structured mesh scheme.

```bash
If you want to use an unstructured mesh scheme, please refer to the EXAMPLE/example7_topography_surface.  
You can generate the mesh by yourself and run the simulation. 
```

## model_elastic_parameters: 

Define velocity and density files for the structured mesh schemes.

```bash
If you want to use an unstructured mesh scheme, you need to define the velocity and density files by 
yourself according to the your mesh schemes. Please refer to the EXAMPLE/example7_topography_surface
```

## openmp: 

An example to show how to use openmp and it is not related to other functions.

## pml: 

Set absorbing boundary condition. 

    I use M-PML developed by myself. You can modify this code to use the PML or C-PML damping profiles, 
    and some corresponding modifications need to be made in the folder: time_evolution/elastic_wave.c. 

## seisfem: 

 Main function to call all other functions to finish the simulation.

## shape: 

 Shape functions for different types of elements. 

## solver: 

Some packages for solving the linear system, including SuperLU_5.2.1 and Paradiso packages. 
But I would like to recommend using the mass-lumped to solve the linear equation. 
Although the mass-lumped technique may cause a long-time simulation unstable problem, 
you can carefully choose the time step to avoid this problem.

```bash
SOLVER_TYPE List:
	I  SOLVER_TYPE    Definition
	-  ------------   ----------
        1  pardiso        Require license and only valid for username: haipeng;
	2  mgmres         Generalized Minimum Residual (GMRES) algorithm, CSR format;
	3  masslump       Mass lump technique.
	
Note: pardiso packages require a license. You need to get your own license and replace pardiso
packages using your own packages and license, which you can get from: https://www.pardiso-project.org
```

## source_receiver: 

Ricker wavelet function. Locations of the source and receiver.      

## sparse_matrix 

Some functions to deal with sparse matrix operations.      	  

## time_evolution 
Solve the elastic wave equation with M-PML, update the wavefiled, and save wavefiled and seismogram files.  


## Note
Please read the README file before you run every example. 
The seisfem code is not that robust and please pay close attention to get things to work well. 
