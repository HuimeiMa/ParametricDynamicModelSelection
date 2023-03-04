# Parametric Dynamic Model Selection
Target equation: the Lorenz 96 System 

udot_{j} = f(t)*( u_{j+1} - u_{j-2} ) * u_{j-1} - u_{j} + \alpha

Auxiliary functions:  
 * BuildMat.m (Restriction of data matrix and velocity vector)
 * legendre.m (The candidate functions)
 * leg2mon.m (Legendre to Monomial Transform)
 * DouglasRachford.m (Optimization routine)
 * Lorenz96Euler.m (Constructs data and velocity vector)
 * SubsetMat.m (Indices used in cyclic permutation and restriction of the data)
 * SupportSet.m (Finds the support set of equation coefficients)
 * dudtFD.m (Velocity vector)  
 * stridge.m (Optimization routine)
 * dudtFD.m (Velocity vector)  
#### Authors: Huimei Ma, Xiaofan Lu, Linan Zhang  
#### Data: September 20, 2022
