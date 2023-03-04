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
sigma = 0.0;% noise level
f = @(t) -(1+cos(t));

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
ctrue = zeros(N,N_t);
ctrue(supp(1),:) = alpha;
ctrue(supp(2),:) = -1; % u_{2}
ctrue(supp(3),:) = f(t); % u_{2} * u_{n}
ctrue(supp(4),:) = -f(t); % u_{n-1} * u_{1,n}

%% Solve for a approximation of ctrue for each t.

Cmon = zeros(N,N_t);
for ii = 1:N_t
    A = Aleg((ii-1)*nb+1:ii*nb,:);
    b1 = V((ii-1)*nb+1:ii*nb,:);
    b = Udot((ii-1)*nb+1:ii*nb,:);
   
    epsilon = 1.2 * norm(b-b1,2); % For testing purposes, in practice must be determined.
    tau = 1; mu = 1/2; MaxIt_1 = 1e5; tol = 1e-6; %Optimization Parameters
    cleg = DouglasRachford(A,b1,sigma,tau,mu,MaxIt_1,tol);

    c = leg2mon(cleg,p,Ind1,Ind20,Ind11,Ind300,Ind210,Ind120,Ind111);
    %Optional Thresholding, See Proposition 3.5
    c = c.*(abs(c)>1e-1);
    Cmon(:, ii) = c;
end

%% Model Identification

tt_start_ind = 2;
tt_end_ind = 30;
tt = t(tt_start_ind:tt_end_ind)';
y1 = Cmon(supp(1),tt_start_ind:tt_end_ind)';
y2 = Cmon(supp(2),tt_start_ind:tt_end_ind)';
y3 = Cmon(supp(3),tt_start_ind:tt_end_ind)';
y4 = Cmon(supp(4),tt_start_ind:tt_end_ind)';

theta = [ones(size(tt)) tt tt.^2 sin(tt) sin(2*tt) cos(tt) cos(2*tt) exp(tt) exp(2*tt)];

MaxIt_2 = 10;
gamma1 = 0.00000001;
lambda1 = 8;
w1 = stridge(theta,y1,gamma1,lambda1,MaxIt_2);
gamma2 = 0.00000001;
lambda2 = 1.1;
w2 = stridge(theta,y2,gamma2,lambda2,MaxIt_2);
gamma3 = 0.000001;
lambda3 =2;
w3 = stridge(theta,y3,gamma3,lambda3,MaxIt_2);
gamma4 = 0.000001;
lambda4 = 2;
w4 = stridge(theta,y4,gamma4,lambda4,MaxIt_2);
w = zeros(36,1);
w(1:9,:) = w1;
w(10:18,:) = w2;
w(19:27,:) = w3;
w(28:36,:) = w4;

w_true = zeros(36,1);
w_true(1,:) = alpha;
w_true(10,:) = -1;
w_true(19,:) = -1;
w_true(24,:) = -1;
w_true(28,:) = 1;
w_true(33,:) = 1;

%% print result

%RMS Error
fprintf('w1 = %f, E(C_K) = %f, E(Gamma_k) = %f\n',w1(1),norm(Cmon(supp(1),tt_start_ind:tt_end_ind)-ctrue(supp(1),tt_start_ind:tt_end_ind))...
    /norm(ctrue(supp(1),tt_start_ind:tt_end_ind),2),norm(w1-w_true(1:9,:))/norm(w_true(1:9,:),2));
fprintf('w2 = %f, E(C_K) = %f, E(Gamma_k) = %f\n',w2(1),norm(Cmon(supp(2),tt_start_ind:tt_end_ind)-ctrue(supp(2),tt_start_ind:tt_end_ind))...
    /norm(ctrue(supp(2),tt_start_ind:tt_end_ind),2),norm(w2-w_true(10:18,:))/norm(w_true(10:18,:),2));
fprintf('w3 = %f,%f, E(C_K) = %f, E(Gamma_k) = %f\n',w3(1),w3(6),norm(Cmon(supp(3),tt_start_ind:tt_end_ind)-ctrue(supp(3),tt_start_ind:tt_end_ind))...
    /norm(ctrue(supp(3),tt_start_ind:tt_end_ind),2),norm(w3-w_true(19:27,:))/norm(w_true(19:27,:),2));
fprintf('w4 = %f, %f, E(C_K) = %f, E(Gamma_k) = %f\n',w4(1),w4(6),norm(Cmon(supp(4),tt_start_ind:tt_end_ind)-ctrue(supp(4),tt_start_ind:tt_end_ind))...
    /norm(ctrue(supp(4),tt_start_ind:tt_end_ind),2),norm(w4-w_true(28:36,:))/norm(w_true(28:36,:),2));
E_C = norm(Cmon(:,tt_start_ind:tt_end_ind) - ctrue(:,tt_start_ind:tt_end_ind))/norm(ctrue(:,tt_start_ind:tt_end_ind))
E_Gamma = norm(w-w_true)/norm(w_true)
SNR = 10*log10(norm(u - mean(u))/norm(tilde_u-u))
