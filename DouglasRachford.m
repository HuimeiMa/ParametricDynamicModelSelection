function x = DouglasRachford(A, b, epsilon, gamma, mu, MaxIt_1, tol)

% ============================================================
% Description:
%   Douglas-Rachford Algorithm
% Inputs:
%   A = m*n matrix
%   b = m*1 vector
%   epsilon = constrain parameter
%   tau = step size
%   mu = convergence parameter
%   MaxIt_1 = maximum number of iterations allowed
%   tol = tolerance
% Outputs:
%   x = min_{x} ||x||_1
%       sub to  ||Ax-b||_2 <= epsilon
% Reference:
%   (1) Hayden Schaeffer, Giang Tran, Rachel Ward, Linan Zhang
%   Extracting structured dynamical systems using sparse optimization with very few samples
%   https://arxiv.org/abs/1805.04158
%   (2) Gabriel Peyre
%   http://www.numerical-tours.com/matlab/sparsity_5_dictionary_learning_denoising/
% Authors: Huimei Ma, Xiaofan Lu, Linan Zhang.
% ============================================================

% Normalize A.
Acnorm = sqrt(sum(A.^2,1)); 
Acnorm = Acnorm(:);
A = normc(A);

% Initialization
n = size(A, 2); 
m = length(b);
x = zeros(n,1); 
tilde_x = zeros(n,1);
tilde_omega = zeros(m,1); 
k = 1; err = tol+1;

pAlower = chol(eye(n) + A'*A,'lower');
mu1 = 1-mu/2;

while k<=MaxIt_1 && err>=tol

    % first step
    [omega_tmp, x_tmp] = rProxG1(tilde_omega, tilde_x);
    [omega_tmp, x_tmp] = rProxG2(omega_tmp, x_tmp);
    tilde_omega = mu1*tilde_omega + mu/2*omega_tmp;
    tilde_x = mu1*tilde_x + mu/2*x_tmp;

    % second step
    [~,x_new] = ProxG1(tilde_omega,tilde_x);

    % update
    err = norm(x_new-x);
    x = x_new;
    k = k+1;

end
x = x./Acnorm;

    function [omega_new, x_new] = ProxG1(omega, x) % ell^1
        x_new = max(abs(x)-gamma, 0) .* sign(x); % ell^1   
        omega_new = b + (omega-b) * min(epsilon/norm(omega-b,2),1);
    end

    function [omega_new,x_new] = ProxG2(omega, x)
        x_new = pAlower'\(pAlower\(x + A'*omega));
        omega_new = A*x_new;
    end

    function [omega_new, x_new] = rProxG1(omega, x)
        [omega_new, x_new] = ProxG1(omega,x);
        omega_new = 2*omega_new-omega; 
        x_new = 2*x_new-x;
    end

    function [omega_new, x_new] = rProxG2(omega, x)
        [omega_new, x_new] = ProxG2(omega, x); 
        omega_new = 2*omega_new-omega; 
        x_new = 2*x_new-x;
    end

end