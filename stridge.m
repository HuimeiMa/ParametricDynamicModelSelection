function c = stridge(Theta, V, gamma, lambda,Maxit)
% ============================================================
% Description:
%   The function solves the linear system Theta * c = V.
% Inputs: 
%   theta = normalizes the columns of D to a length of 1.
%   y = velocity vector
%   gamma = parameters of Cholesky's decomposition
%   lambda = parameters of sparsification
%   theta_cnorm = the norm of each column.
% Outputs:
%   c = parametric coefficient
%
% Reference: Data-Driven Discovery of PDEs
%            https://www.science.org/doi/10.1126/sciadv.1602614
%
% Authors: Samuel Rudy, Alessandro Alla, Steven L. Brunton, J. Nathan Kutz
% ============================================================

% Column-normalize Theta.
theta_cnorm_factor = sqrt(sum(Theta.^2,1));
theta_cnorm_factor = theta_cnorm_factor(:);
Theta_cnorm = normc(Theta);

% STRidge
c = inverse(Theta_cnorm, V);
[~,n] = size(c);
for k = 1:Maxit
    smallinds = (abs(c) < lambda);
    c(smallinds) = 0;
    biginds = ~smallinds;
    for ind = 1:n
        c(biginds, ind) = inverse(Theta_cnorm(:,biginds), V(:,ind));
    end
end
c = c./theta_cnorm_factor;

    function x = inverse(A, b)
        A_chol = chol(gamma*eye(size(A,2)) + A'*A,'lower');
        x = A_chol'\(A_chol\(A'*b));
    end
end