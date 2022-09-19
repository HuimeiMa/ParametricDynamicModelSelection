function [U, V, Udot] = BuildMat(u, ut, udot, t, indr, indc)

% ============================================================
% Description:
%   The function outputs restriction of data matrix and velocity vector
%
% Inputs:
%   u = the state of a nonlinear system, i.e. given data
%   ut = 1th time derivative of u
%   indt = time index
%   indr = row index
%   indc = column index (restriction of the data)
%   Nx = number of sampling space
%   nt = The length of time index
%   nx = The length of column index (The length of row index )
%
% Outputs:
%   U = restriction of data matrix
%   V = velocity vector -- computed velocity
%
% Example: Nx=5, nt=2 nb=3, ii=1, 
%           data matrix U = | u1(t_0) u2(t_0) u3(t_0) u4(t_0) u5(t_0) |
%                           | u1(t_1) u2(t_1) u3(t_1) u4(t_1) u5(t_1) |
%           velocity vector V = | v1(t_0) v2(t_0) v3(t_0) v4(t_0) v5(t_0) |
%                               | v1(t_1) v2(t_1) v3(t_1) v4(t_1) v5(t_1) |
%   cyclic permutation of the data =>
%               | u1(t_0) u2(t_0) u3(t_0) u4(t_0) u5(t_0) |
%               | u2(t_0) u3(t_0) u4(t_0) u5(t_0) u1(t_0) |
%               | u3(t_0) u4(t_0) u5(t_0) u1(t_0) u2(t_0) |
%               | u4(t_0) u5(t_0) u1(t_0) u2(t_0) u3(t_0) |
%               | u5(t_0) u1(t_0) u2(t_0) u3(t_0) u4(t_0) |
%           U = | u1(t_1) u2(t_1) u3(t_1) u4(t_1) u5(t_1) |
%               | u2(t_1) u3(t_1) u4(t_1) u5(t_1) u1(t_1) |
%               | u3(t_1) u4(t_1) u5(t_1) u1(t_1) u2(t_1) |
%               | u4(t_1) u5(t_1) u1(t_1) u2(t_1) u3(t_1) |
%               | u5(t_1) u1(t_1) u2(t_1) u3(t_1) u4(t_1) |
%   restriction of the data =>
%               | u1(t_0) u2(t_0) u5(t_0) |
%               | u2(t_0) u3(t_0) u1(t_0) |
%               | u3(t_0) u4(t_0) u2(t_0) |
%           U = | u1(t_1) u2(t_1) u5(t_1) |
%               | u2(t_1) u3(t_1) u1(t_1) |
%               | u3(t_1) u4(t_1) u2(t_1) |
%   corresponding velocity vector =>
%               | v1(t_0) |
%               | v2(t_0) |
%           V = | v3(t_0) |
%               | v1(t_1) |
%               | v2(t_1) |
%               | v3(t_1) |
% Reference: Extracting Structured Dynamical Systems Using Sparse Optimization
%            With Very Few Samples
%            https://epubs.siam.org/doi/10.1137/18M1194730
%
% Authors: Hayden Schaeffer, Giang Tran, Rachel Ward, Linan Zhang
% ============================================================
Nx = size(u,2);
nx = length(indc);
nt = length(t);
U = zeros(nt*nx, nx); % data matrix
V = zeros(nt*nx, 1); % velocity vector

for k = 1:nt
    % Perform cyclic permutation on U.
    v = zeros(nx, Nx); % temporal variable
    v(1,:) = circshift(u(k,:),[0 -indr(1)+1]);
    for j=2:nx
        v(j,:) = circshift(v(j-1,:),[0 -1]);
    end

    % Restrict the data onto the block.
    U((k-1)*nx+1 : (k-1+1)*nx, :) = v(:,indc);
    % Construct  V accordingly
    V((k-1)*nx+1 : (k-1+1)*nx, :) = ut(k, indr)';
    Udot((k-1)*nx+1 : (k-1+1)*nx, :) = udot(k, indr)';
end
end