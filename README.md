# ParametricDynamicModelSelection
Target equation: the lorenz 96 system 
$$\dot{u}_{i} =f(t) ( u_{i+1} - u_{i-2} )  u_{i-1} - u_{i} + \alpha$$  
Auxiliary functions:
 * Lorenz96Euler.m(Constructs data and velocity vector)
 * SubsetMat.m(Indices used in cyclic permutation and restriction of the data)
 * BuildMat.m(Restriction of data matrix and velocity vector)
 * Dictionary.m(Dictioary matrix)
 * SupportSet.m(Finds the support set of equation coefficients)
 * DouglasRachford.m(Optimization Routine)
 * stridge.m(Optimization Routine)
 * dudtFD.m(Velocity vector)  
#### Authors: Huimei Ma, Xiaofan Lu, Linan Zhang  
#### Data: September 20, 2022
