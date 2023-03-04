function [u,udot,udot1] = Lorenz96Euler(N,t,f,alpha,sigma)

% ============================================================
% Description:
%   The function constructs data and velocity vector
% Inputs:
%   n = number of grid points in each direction
%   t = temporal interval
%   f,F = Parameters of Lorenz equation
% Outputs:
%   u = the state of the ODE system, size(u) = Nt-by-n
%   udot = du/dt, size(udot) = Nt-by-n
%   udot1 = numerical derivative
%
% Authors: Huimei Ma, Xiaofan Lu, Linan Zhang.
% ============================================================
u0 = 2*rand(N,1)-1;
Nt = length(t);
delta_t = t(2) - t(1);
% Refine the interval.
% t1,t2 -> t1,t11,t12,...,t19,t2
np = 1000;
ddt = delta_t/np;

% Initialize u and udot.
d = N;
nt = round((t(end)-t(1))/ddt)+1;
tt = linspace(t(1),t(end),nt);

% Update u and udot.
u = zeros(Nt,d);
udot = zeros(Nt-1,d);
u2 = u0;
udot2 = dudt(u2,t(1),alpha);
u(1,:) = u2; 
udot(1,:) = udot2;
% update u(2:nt) and udot(1:nt-1)
for ii=2:nt
    % ii = r + np * jj
    r = mod(ii,np); jj = floor(ii/np);
    u2 = u2 + udot2*ddt;
    udot2 = dudt(u2, tt(ii), alpha);
    if r==1
        u(jj+1,:) = u2;
        udot(jj+1,:) = udot2;
    end
end

% Compute numerical derivative.
udot1 = dudtFD(u,delta_t);
    function udot = dudt(u,t_i,F)
        udot = f(t_i)*(circshift(u,-1)-circshift(u,2)) ...
            .* circshift(u,1) - u + F;
    end


end
