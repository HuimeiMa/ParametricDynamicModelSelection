function x = DouglasRachford(A, b, sigma, tau, mu, MaxIt, tol)

% ============================================================
% Description:
%   Douglas-Rachford Algorithm
% Inputs:
%   A = m*n matrix
%   b = m*1 vector
%   sigma = constrain parameter
%   tau = step size
%   mu = convergence parameter
%   MaxIt = maximum number of iterations allowed
%   tol = tolerance
% Outputs:
%   u = min_{x} ||x||_1
%       sub to  ||Ax-b||_2 <= sigma
% Reference:
%   (1) Hayden Schaeffer, Giang Tran, Rachel Ward, Linan Zhang
%   Extracting structured dynamical systems using sparse optimization with very few samples
%   https://arxiv.org/abs/1805.04158
%   (2) Gabriel Peyre
%   http://www.numerical-tours.com/matlab/sparsity_5_dictionary_learning_denoising/
% Authors: Huimei Ma, Xiaofan Lu, Linan Zhang.
% ============================================================

% Normalize A.
Acnorm = sqrt(sum(A.^2,1)); Acnorm = Acnorm(:);
A = normc(A);

% Initialization
N = size(A, 2); M = length(b);
x = zeros(N,1); x1 = zeros(N,1); u1 = zeros(M,1);
T = 1; error = tol+1;

pAlower = chol(eye(N) + A'*A,'lower');
mu1 = 1-mu;

while T<=MaxIt && error>=tol

    % first step
    [p, q] = rProxF(x1, u1);
    [p, q] = rProxG(p, q);
    x1 = mu1*x1 + mu*p;
    u1 = mu1*u1 + mu*q;

    % second step
    [xnew,~] = ProxF(x1,u1);

    % update
    error = norm(x-xnew);
    x = xnew;
    T = T+1;

end
x = x./Acnorm;

    function [x1, u1] = ProxF(x, u) % ell^1
        x1 = max(abs(x)-tau, 0) .* sign(x); % ell^1
        u1 = b + (u-b) * min(sigma/norm(u-b,2),1);
    end

    function [x1,u1] = ProxG(x, u)
        x1 = pAlower'\(pAlower\(x + A'*u));
        u1 = A*x1;
    end

    function [x1, u1] = rProxF(x, u)
        [x1, u1] = ProxF(x,u); x1 = 2*x1-x; u1 = 2*u1-u;
    end

    function [x1, u1] = rProxG(x, u)
        [x1, u1] = ProxG(x, u); x1 = 2*x1-x; u1 = 2*u1-u;
    end

end