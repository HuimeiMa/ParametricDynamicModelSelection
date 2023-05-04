clc; clear all; close all;
set(0,'defaultaxesfontsize',13);
% =================================================================
% Target equation: the Lorenz 96 equation
%       udot_{j} = f(t)*( u_{j+1} - u_{j-2} ) * u_{j-1} - u_{j} + F
%       f = @(t) -(1+cos(t));
% =================================================================

%% Create the entire data.

% parameters
N = 128;
alpha = 8;
delta_t = 0.002;
t = 0:delta_t:(100)*delta_t;
N_t = length(t);
sigma = 0e-06;% noise level
f = @(t) 3*exp(t);%f(t)=-(1+cos(t)) or f(t)=3*exp(t)

%Constructed data and velocity vector
[u,udot,~] = Lorenz96Euler(N,t,f,alpha,sigma);
tilde_u = u + sigma*randn(size(u));
udot1 = dudtFD(tilde_u,delta_t);
%% Downsample the data.

%Outputs indices used in cyclic permutation and restriction of the data
ii = 1; % the index in the entire data for u1 of the block 
nb = 55;% size of the block
[indr, indc] = SubsetMat(ii, nb, N); 

%Outputs restriction of data matrix and velocity vector
[U, V, Udot] = BuildMat(tilde_u, udot1, udot, t, indr, indc);
%% Create the dictionary.

p = 2; % degree of the basis element (p=2 or p=3)
r = 5; % localization of the dictionary (radius of the restricted subset)
[Amon,Aleg,L,Ind1,Ind20,Ind11,Ind300,Ind210,Ind120,Ind111] = legendre(U,p,r,indc);
N = size(Aleg,2);   

% support set
supp = SupportSet(L, indc, 'lorenz96');

% Exact coefficients
Ctrue = zeros(N,N_t);
Ctrue(supp(1),:) = alpha;
Ctrue(supp(2),:) = -1; % u_{2}
Ctrue(supp(3),:) = f(t); % u_{2} * u_{n}
Ctrue(supp(4),:) = -f(t); % u_{n-1} * u_{1,n}

%% Solve for a approximation of ctrue for each t.

Cmon = zeros(N,N_t);
for ii = 1:N_t
    Phi = Aleg((ii-1)*nb+1:ii*nb,:);
    b1 = V((ii-1)*nb+1:ii*nb,:);
    b = Udot((ii-1)*nb+1:ii*nb,:);
   
    epsilon = 1.2 * norm(b-b1,2); % For testing purposes, in practice must be determined.
    gamma = 1; mu = 1/2; MaxIt_1 = 1e5; tol = 1e-6; %Optimization Parameters
    cleg = DouglasRachford(Phi,b1,sigma,gamma,mu,MaxIt_1,tol);

    C = leg2mon(cleg,p,Ind1,Ind20,Ind11,Ind300,Ind210,Ind120,Ind111);
    %Optional Thresholding, See Proposition 3.5
    C = C.*(abs(C)>1e-1);
    Cmon(:, ii) = C;
end

%% Model Identification
tt_start_ind = 2;
tt_end_ind = 30;
tt = t(tt_start_ind:tt_end_ind)';

Psi= [ones(size(tt)) tt tt.^2 sin(tt) sin(2*tt) cos(tt) cos(2*tt) exp(tt) exp(2*tt)];
%If the reader wants to learn the lower-order Taylor approximation of the governing equations, 
% the library Psi is available:
%Psi = [ones(size(tt)) tt tt.^2 tt.^3];
[w,w_true] = ModelIdentification(tt_start_ind,tt_end_ind,alpha,Psi,Cmon,Ctrue,supp,'exp(t)');
%ft_name = 'cos(t) or 'exp(t)'


%% print result

E_C = norm(Cmon(:,tt_start_ind:tt_end_ind) - Ctrue(:,tt_start_ind:tt_end_ind))/norm(Ctrue(:,tt_start_ind:tt_end_ind))
E_Gamma = norm(w-w_true)/norm(w_true)
SNR = 10*log10(norm(u - mean(u))/norm(tilde_u-u))
