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
[Phi,L] = Dictionary(U,p,r,indc);
N = size(Phi,2);    

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
    A = Phi((ii-1)*nb+1:ii*nb,:);
    b1 = V((ii-1)*nb+1:ii*nb,:);
    b = Udot((ii-1)*nb+1:ii*nb,:);
   
    epsilon = 1.2 * norm(b-b1,2); % For testing purposes, in practice must be determined.
    tau = 1; mu = 1/2; MaxIt_1 = 1e5; tol = 1e-6; %Optimization Parameters
    c = DouglasRachford(A,b1,epsilon,tau,mu,MaxIt_1,tol);
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

%% print result

%RMS Error
fprintf('w1 = %f, RMS error = %f\n',w1(1),norm(Cmon(supp(1),tt_start_ind:tt_end_ind)-ctrue(supp(1),tt_start_ind:tt_end_ind))...
    /norm(ctrue(supp(1),tt_start_ind:tt_end_ind),2));
fprintf('w2 = %f, RMS error = %f\n',w2(1),norm(Cmon(supp(2),tt_start_ind:tt_end_ind)-ctrue(supp(2),tt_start_ind:tt_end_ind))...
    /norm(ctrue(supp(2),tt_start_ind:tt_end_ind),2));
fprintf('w3 = %f,%f, RMS error = %f\n',w3(1),w3(6),norm(Cmon(supp(3),tt_start_ind:tt_end_ind)-ctrue(supp(3),tt_start_ind:tt_end_ind))...
    /norm(ctrue(supp(3),tt_start_ind:tt_end_ind),2));
fprintf('w4 = %f, %f, RMS error = %f\n',w4(1),w4(6),norm(Cmon(supp(4),tt_start_ind:tt_end_ind)-ctrue(supp(4),tt_start_ind:tt_end_ind))...
    /norm(ctrue(supp(4),tt_start_ind:tt_end_ind),2));


% cosXi = figure(2)
% subplot(2,2,1)
% plot(t(2:end-1),ctrue(supp(1),2:end-1),'color','black','LineWidth',1.5)
% hold on 
% plot(t(2:end-1),Cmon(supp(1),2:end-1),'*')
% xlabel('t')
% xline(0.004,'--');
% xline(0.06,'--');
% legend('True coefficients','Approximation','Location','southeast')
% title('Coefficients of the term F as a function of t')
% set(gca,'Xlim',[0,0.2],'YLim',[7.95,8.05]);
% 
% subplot(2,2,2)
% plot(t(2:end-1),ctrue(supp(2),2:end-1),'color','black','LineWidth',1.5)
% hold on 
% plot(t(2:end-1),Cmon(supp(2),2:end-1),'*')
% xlabel('t')
% xline(0.004,'--');
% xline(0.06,'--');
% legend('True coefficients','Approximation','Location','southeast')
% title('Coefficients of the term u_j as a function of t')
% set(gca,'YLim',[-1.05,-0.95]);
% 
% subplot(2,2,3)
% plot(t(2:end-1),ctrue(supp(3),2:end-1),'color','black','LineWidth',1.5)
% hold on 
% plot(t(2:end-1),Cmon(supp(3),2:end-1),'*')
% xlabel('t')
% xline(0.004,'--');
% xline(0.06,'--');
% legend('True coefficients','Approximation','Location','southeast')
% title('Coefficients of the term u_{j-2}u_{j-1} as a function of t')
% set(gca,'YLim',[-2.01,-1.98]);
% 
% subplot(2,2,4)
% plot(t(2:end-1),ctrue(supp(4),2:end-1),'color','black','LineWidth',1.5)
% hold on 
% plot(t(2:end-1),Cmon(supp(4),2:end-1),'*')
% xlabel('t')
% xline(0.004,'--');
% xline(0.06,'--');
% legend('True coefficients','Approximation','Location','southeast')
% title('Coefficients of the term u_{j+1}u_{j-1} as a function of t')
% set(gca,'YLim',[1.97,2.00]);
% 
% fig = gcf;
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 12 6];
% print(cosXi,'-dpng','-r0');