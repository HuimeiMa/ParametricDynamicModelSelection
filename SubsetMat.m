function [indr, indc] = SubsetMat(start_index, length_block, Nx)

% ============================================================
% Description:
%   The function outputs indices used in cyclic permutation 
%   and restriction of the data
%
% Inputs:
%   start_index = index of the first point of the block in the entire data
%   length_block = size of the block 
%             Requirements for length_block: length_block odd and
%             length_block+2 <= Nx, and length_block + start_index < Nx
%   Nx = Nx = number of sampling space
% Outputs:
%   indr = row index
%   indc = column index (restriction of the data)
%
% Example: Nx=5, start_index = 1, length_block = 3
%          indr = [1, 2, 3]
%          indc = [1, 2, 5]
% Reference: Extracting Structured Dynamical Systems Using Sparse Optimization
%            With Very Few Samples
%            https://epubs.siam.org/doi/10.1137/18M1194730
%
% Authors: Hayden Schaeffer, Giang Tran, Rachel Ward, Linan Zhang
% ============================================================

r = (length_block - 1)/2;

% indices used in cyclic permutation and restriction of the data
indc = 1:Nx; 
indc = [ indc(1:r+1) indc(Nx-r+1:Nx)];
indc = sort(indc); % column indices

indr = 1:Nx; 
indr = indr(start_index:start_index+length_block-1); % row indices

end