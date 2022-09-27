# Parametric Dynamic Model Selection
Target equation: the Lorenz 96 System 

udot_{i} = f(t)*( u_{i+1} - u_{i-2} ) * u_{i-1} - u_{i} + \alpha

Auxiliary functions:  
 * BuildMat.m (Restriction of data matrix and velocity vector)
 * Dictionary.m (The candidate functions)
 * DouglasRachford.m (Optimization Routine)
 * Lorenz96Euler.m (Constructs data and velocity vector)
 * SubsetMat.m (Indices used in cyclic permutation and restriction of the data)
 * SupportSet.m (Finds the support set of equation coefficients)
 * dudtFD.m (Velocity vector)  
 * stridge.m (Optimization Routine)
 * dudtFD.m (Velocity vector)  
#### Authors: Huimei Ma, Xiaofan Lu, Linan Zhang  
#### Data: September 20, 2022
