function x = STRidge(A, b, lambda_2, lambda_1, Maxit_2)
% ============================================================
% Description:
%   The function solves the linear system Theta * c = V.
% Inputs: 
%   A = normalizes the columns of D to a length of 1.
%   b = velocity vector
%   lambda_2 = parameters of Cholesky's decomposition
%   lambda_1 = parameters of sparsification
%   A_cnorm = the norm of each column.
% Outputs:
%   x = parametric coefficient
%
% Reference: Data-Driven Discovery of PDEs
%            https://www.science.org/doi/10.1126/sciadv.1602614
%
% Authors: Samuel Rudy, Alessandro Alla, Steven L. Brunton, J. Nathan Kutz
% ============================================================

% Column-normalize Theta.
A_cnorm_factor = sqrt(sum(A.^2,1));
A_cnorm_factor = A_cnorm_factor(:);
A_cnorm = normc(A);

% STRidge
x = inverse(A_cnorm, b);
[~,n] = size(x);
for k = 1:Maxit_2
    smallinds = (abs(x) < lambda_1);
    x(smallinds) = 0;
    biginds = ~smallinds;
    for ind = 1:n
        x(biginds, ind) = inverse(A_cnorm(:,biginds), b(:,ind));
    end
end
x = x./A_cnorm_factor;

    function x = inverse(A, b)
        A_chol = chol(lambda_2*eye(size(A,2)) + A'*A,'lower');
        x = A_chol'\(A_chol\(A'*b));
    end
end