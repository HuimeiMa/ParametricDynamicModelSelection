clc; clear all; close all;
set(0,'defaultaxesfontsize',11);
% =================================================================
% Target equation: the Lorenz 96 equation
%       udot_{i} = f(t)*( u_{i+1} - u_{i-2} ) * u_{i-1} - u_{i} + F
%       f = @(t) -(1+cos(t));
% =================================================================

% parameters
n = 128;
F = 8;
dt = 0.002;
t = 0:dt:100*dt;
Nt = length(t);
var = 0.0001;% noise level
f = @(t) -(1+cos(t));

%Constructed data and velocity vector
[u,udot,udot1] = Lorenz96Euler(n,t,f,F);

%% Downsample the data.

%Outputs indices used in cyclic permutation and restriction of the data
ii = 1; % the index in the entire data for u1 of the block 
nb = 55;% size of the block
[indr, indc] = SubsetMat(ii, nb, n); 

%Outputs restriction of data matrix and velocity vector
[U, V, Udot] = BuildMat(u, udot1, udot, t, indr, indc);
%% Create the dictionary.

p = 2; % degree of the basis element (p=2 or p=3)
r = 5; % localization of the dictionary (radius of the restricted subset)
[D,L] = Dictionary(U,p,r,indc);
N = size(D,2);    

% support set
supp = SupportSet(L, indc, 'lorenz96');

% Exact coefficients
ctrue = zeros(N,Nt);
ctrue(supp(1),:) = F;
ctrue(supp(2),:) = -1; % u_{2}
ctrue(supp(3),:) = f(t); % u_{2} * u_{n}
ctrue(supp(4),:) = -f(t); % u_{n-1} * u_{1,n}

%% Solve for a approximation of ctrue for each t.

Cmon = zeros(N,Nt);
for ii = 1: Nt
    A = D((ii-1)*nb+1:ii*nb,:);
    b1 = V((ii-1)*nb+1:ii*nb,:);
    b = Udot((ii-1)*nb+1:ii*nb,:);
   
    sigma = 1.01 * norm(b-b1,2); % For testing purposes, in practice must be determined.
    tau = 1; mu = 1/2; MaxIt = 1e5; tol = 1e-6; %Optimization Parameters
    c = DouglasRachford(A,b1,sigma,tau,mu,MaxIt,tol);
    Cmon(:, ii) = c;
end

%% Model Identification

tt_start_ind = 2;
tt_end_ind = 30;
tt = t(tt_start_ind:tt_end_ind);
y1 = Cmon(supp(1),tt_start_ind:tt_end_ind)';
y2 = Cmon(supp(2),tt_start_ind:tt_end_ind)';
y3 = Cmon(supp(3),tt_start_ind:tt_end_ind)';
y4 = Cmon(supp(4),tt_start_ind:tt_end_ind)';

theta = [];
for i = 1:4
    theta = [theta tt'.^(i-1)]; 
end

MaxIt2 = 10;
gamma1 = 0.0000001;
lambda1 = 10;
w1 = stridge(theta,y1,gamma1,lambda1,MaxIt2);
gamma2 = 0.0000001;
lambda2 = 2;
w2 = stridge(theta,y2,gamma2,lambda2,MaxIt2);
gamma3 = 0.000001;
lambda3 = 0.001;
w3 = stridge(theta,y3,gamma3,lambda3,MaxIt2);
gamma4 = 0.000001;
lambda4 = 0.001;
w4 = stridge(theta,y4,gamma4,lambda4,MaxIt2);

%% print result

%RMS Error
fprintf('w1 = %f, RMS error = %f\n',w1(1),norm(Cmon(supp(1),tt_start_ind:tt_end_ind)-ctrue(supp(1),tt_start_ind:tt_end_ind))...
    /norm(ctrue(supp(1),tt_start_ind:tt_end_ind),2));
fprintf('w2 = %f, RMS error = %f\n',w2(1),norm(Cmon(supp(2),tt_start_ind:tt_end_ind)-ctrue(supp(2),tt_start_ind:tt_end_ind))...
    /norm(ctrue(supp(2),tt_start_ind:tt_end_ind),2));
fprintf('w3 = %f,%f, RMS error = %f\n',w3(1),w3(3),norm(Cmon(supp(3),tt_start_ind:tt_end_ind)-ctrue(supp(3),tt_start_ind:tt_end_ind))...
    /norm(ctrue(supp(3),tt_start_ind:tt_end_ind),2));
fprintf('w4 = %f, %f, RMS error = %f\n',w4(1),w4(3),norm(Cmon(supp(4),tt_start_ind:tt_end_ind)-ctrue(supp(4),tt_start_ind:tt_end_ind))...
    /norm(ctrue(supp(4),tt_start_ind:tt_end_ind),2));

% figure(1);
% [t1,xx] = meshgrid(t,1:128);
% surf(t1,xx,u');
% shading interp
% xlabel('t');%x
% ylabel('j');%y
% zlabel('u');%z
% set(gca,'FontSize',11,'LineWidth',1);

% figure(3)
% plot(tt,ctrue(supp(1),tt_start_ind:tt_end_ind),'LineWidth',2)
% hold on 
% plot(tt,Cmon(supp(1),tt_start_ind:tt_end_ind),'d')
% legend('True coefficients','approximation')
% xlabel('t')
% ylabel('C(t)')
% set(gca,'YLim',[7.95,8.05]);
% 
% figure(3)
% plot(tt,ctrue(supp(2),tt_start_ind:tt_end_ind),'LineWidth',2)
% hold on 
% plot(tt,Cmon(supp(2),tt_start_ind:tt_end_ind),'d')
% legend('True coefficients','approximation')
% xlabel('t')
% ylabel('C(t)')
% set(gca,'YLim',[-1.05,-0.95]);
% 
% figure(4)
% plot(tt,ctrue(supp(3),tt_start_ind:tt_end_ind),'LineWidth',2)
% hold on 
% plot(tt,Cmon(supp(3),tt_start_ind:tt_end_ind),'d')
% legend('True coefficients','approximation')
% xlabel('t')
% ylabel('C(t)')
% 
% figure(5)
% plot(tt,ctrue(supp(4),tt_start_ind:tt_end_ind),'LineWidth',2)
% hold on 
% plot(tt,Cmon(supp(4),tt_start_ind:tt_end_ind),'d')
% legend('True coefficients','approximation')
% xlabel('t')
% ylabel('C(t)')