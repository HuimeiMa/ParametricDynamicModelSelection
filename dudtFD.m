function udot = dudtFD(u,dt)
% ============================================================
% Description:
%   The function outputs velocity vector
%
% Inputs:
%   u = the state of a nonlinear system, i.e. given data
%   dt = time step
%
% Outputs:
%   udot = 1th time derivative of u
% Authors: Huimei Ma, Xiaofan Lu, Linan Zhang.
% ============================================================

udot = zeros(size(u));

udot(1,:) = ( u(2,:) - u(1,:) )/(dt);
if size(u,1)>2
    udot(2:end-1,:) = ( u(3:end,:) - u(1:end-2,:) )/(2*dt);
end
udot(end,:) = ( u(end,:) - u(end-1,:) )/(dt);
end