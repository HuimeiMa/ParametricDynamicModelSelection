function [Dmon,L] = Dictionary(u,p,r,indc)

% ============================================================
% Description:
%   The function outputs the Dictionary Matrix.
%
% Inputs:
%   u = the state of a nonlinear system, i.e. given data
%   p = maximum degree of polynomials in the dictionary (p<=3)
%   r = number of left/right neighbors
%   indc = column index (restriction of the data)
%
% Outputs:
%   Dmon = dictionary of monomials of degree at most p
%   L = legend of D
%
% Example: u=[u1,u2,u3,u8,u9], p=2, r=1, indc=[1,2,3,8,9]
%   Dmon = [1, u1, u2, u9, u1^2, u1u2, u1u9, u2^2, u2u9, u9^2]
%
% Authors: Hayden Schaeffer, Giang Tran, Rachel Ward, Linan Zhang
% Date: May 9, 2018
% ============================================================

% Compute the size of u.
% m = number of time grid
% n = number of states.
[m,n] = size(u);
% List the name of u1.
l1 = cell(1,n+1);
for ii=1:n
    l1{ii+1} = strcat('u',num2str(indc(ii)));
end
% Augment u to [1;u] and down-sample the columns.
if n>= 2*r+1
    umon = [ones(m,1) u(:,1:r+1) u(:,n-r+1:n)];
    uleg = [ones(m,1) u(:,1:r+1)*sqrt(3) u(:,n-r+1:n)*sqrt(3)];
    l1 = l1([1:r+2,n-r+2:n+1]);
else
    umon = [ones(m,1) u];
    uleg = [ones(m,1) u*sqrt(3)];
end

% Redefine n.
n = min(2*r+1,n);

% Find all polynomials of degree at most 3.
C = 1:(n+1);
for ii=1:p-1
    C = combvec(C,1:(n+1));
end
% Sort each column of C in ascending order.
C = sort(C,1);
% Remove duplicate rows in C' (i.e. duplicate columns in C).
C = unique(C','rows');
% number of rows of C = number of columns of D
% Note: number of columns of C = 3
nD = size(C,1);
% Construct D.
Dmon = ones(m,nD);
L = cell(nD,1);

%   P2: 3*ui*uj for i~=j
%   P3: sqrt(27)*ui*uj*uk
for ii=1:nD
    % Multiply the corresponding columns to make a polynomial of degree at most p.
    for jj = 1:p
        Dmon(:,ii) = Dmon(:,ii) .* umon(:,C(ii,jj));
        L{ii} = strcat(L{ii},l1{C(ii,jj)});
    end
end
L{1} = num2str(1);
end