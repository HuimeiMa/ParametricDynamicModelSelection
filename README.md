# Parametric Dynamic Model Selection
Target equation: the Lorenz 96 System 

udot_{i} = f(t)*( u_{i+1} - u_{i-2} ) * u_{i-1} - u_{i} + \alpha

Auxiliary functions:   
 * Dictionary.m (The Candidate Functions)
 * SubsetMat.m (Indices used in cyclic permutation and restriction of the data)
 * BuildMat.m (Restriction of data matrix and velocity vector)
 * SupportSet.m (Finds the support set of equation coefficients)
 * DouglasRachford.m (Optimization Routine)
 * stridge.m (Optimization Routine)
 * dudtFD.m (Velocity vector)  
 * Lorenz96Euler.m (constructs data and velocity vector)
#### Authors: Huimei Ma, Xiaofan Lu, Linan Zhang  
#### Data: September 20, 2022
