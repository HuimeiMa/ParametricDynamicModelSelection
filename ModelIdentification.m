function [w,w_true] = ModelIdentification(tt_start_ind,tt_end_ind,alpha,Psi,Cmon,Ctrue,supp,ft_name)
% ============================================================
% Description:
%   The function obtains the final expression of the governing equation.
%
% Inputs:
%   tt_start_ind =  Index of the start time point
%   tt_end_ind = Index of the end point in time
%   alpha = Coefficient of the Lorenz 96 System
%   Cmon =  The coefficient matrix for each time stamp
%   Ctrue = The true coefficient matrix for each time stamp
%   supp = The support set
%   ft_name = The symbol of a time-varying parameter
% Outputs:
%   w = The coefficient of a time-varying parameter
%   w_true = The true coefficient of a time-varying parameter
%
% Authors: Huimei Ma, Xiaofan Lu, Linan Zhang.
% ============================================================

y1 = Cmon(supp(1),tt_start_ind:tt_end_ind)';
y2 = Cmon(supp(2),tt_start_ind:tt_end_ind)';
y3 = Cmon(supp(3),tt_start_ind:tt_end_ind)';
y4 = Cmon(supp(4),tt_start_ind:tt_end_ind)';

[~,m] = size(Psi);

if m>=9
    if isequal('cos(t)', ft_name)
        MaxIt_2 = 10;
        gamma1 = 1e-08;
        lambda1 = 8;
        w1 = STRidge(Psi,y1,gamma1,lambda1,MaxIt_2);
        gamma2 = 1e-08;
        lambda2 = 1.1;
        w2 = STRidge(Psi,y2,gamma2,lambda2,MaxIt_2);
        gamma3 = 1e-06;
        lambda3 =2;
        w3 = STRidge(Psi,y3,gamma3,lambda3,MaxIt_2);
        gamma4 = 2e-06;
        lambda4 = 2.2;
        w4 = STRidge(Psi,y4,gamma4,lambda4,MaxIt_2);
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
        %RMS Error
        fprintf('w1 = %f, E(C_K) = %f, E(Gamma_k) = %f\n',w1(1),norm(Cmon(supp(1),tt_start_ind:tt_end_ind)-Ctrue(supp(1),tt_start_ind:tt_end_ind))...
            /norm(Ctrue(supp(1),tt_start_ind:tt_end_ind),2),norm(w1-w_true(1:9,:))/norm(w_true(1:9,:),2));
        fprintf('w2 = %f, E(C_K) = %f, E(Gamma_k) = %f\n',w2(1),norm(Cmon(supp(2),tt_start_ind:tt_end_ind)-Ctrue(supp(2),tt_start_ind:tt_end_ind))...
            /norm(Ctrue(supp(2),tt_start_ind:tt_end_ind),2),norm(w2-w_true(10:18,:))/norm(w_true(10:18,:),2));
        fprintf('w3 = %f,%f, E(C_K) = %f, E(Gamma_k) = %f\n',w3(1),w3(6),norm(Cmon(supp(3),tt_start_ind:tt_end_ind)-Ctrue(supp(3),tt_start_ind:tt_end_ind))...
            /norm(Ctrue(supp(3),tt_start_ind:tt_end_ind),2),norm(w3-w_true(19:27,:))/norm(w_true(19:27,:),2));
        fprintf('w4 = %f, %f, E(C_K) = %f, E(Gamma_k) = %f\n',w4(1),w4(6),norm(Cmon(supp(4),tt_start_ind:tt_end_ind)-Ctrue(supp(4),tt_start_ind:tt_end_ind))...
            /norm(Ctrue(supp(4),tt_start_ind:tt_end_ind),2),norm(w4-w_true(28:36,:))/norm(w_true(28:36,:),2));

    elseif isequal('exp(t)', ft_name)
        MaxIt_2 =10;
        gamma1 = 1e-08;
        lambda1 = 8;
        w1 = STRidge(Psi,y1,gamma1,lambda1,MaxIt_2);
        gamma2 = 1e-08;
        lambda2 = 1;
        w2 = STRidge(Psi,y2,gamma2,lambda2,MaxIt_2);
        gamma3 = 1e-08;
        lambda3 = 3;
        w3 = STRidge(Psi,y3,gamma3,lambda3,MaxIt_2);
        gamma4 = 1e-08;
        lambda4 = 3;
        w4 = STRidge(Psi,y4,gamma4,lambda4,MaxIt_2);

        w = zeros(36,1);
        w(1:9,:) = w1;
        w(10:18,:) = w2;
        w(19:27,:) = w3;
        w(28:36,:) = w4;

        w_true = zeros(36,1);
        w_true(1,:) = alpha;
        w_true(10,:) = -1;
        w_true(26,:) = 3;
        w_true(35,:) = -3;
        %RMS Error
        fprintf('w1 = %f, E(C_K) = %f, E(Gamma_k) = %f\n',w1(1),norm(Cmon(supp(1),tt_start_ind:tt_end_ind)-Ctrue(supp(1),tt_start_ind:tt_end_ind))...
            /norm(Ctrue(supp(1),tt_start_ind:tt_end_ind),2),norm(w1-w_true(1:9,:))/norm(w_true(1:9,:),2));
        fprintf('w2 = %f, E(C_K) = %f, E(Gamma_k) = %f\n',w2(1),norm(Cmon(supp(2),tt_start_ind:tt_end_ind)-Ctrue(supp(2),tt_start_ind:tt_end_ind))...
            /norm(Ctrue(supp(2),tt_start_ind:tt_end_ind),2),norm(w2-w_true(10:18,:))/norm(w_true(10:18,:),2));
        fprintf('w3 = %f, E(C_K) = %f, E(Gamma_k) = %f\n',w3(8),norm(Cmon(supp(3),tt_start_ind:tt_end_ind)-Ctrue(supp(3),tt_start_ind:tt_end_ind))...
            /norm(Ctrue(supp(3),tt_start_ind:tt_end_ind),2),norm(w3-w_true(19:27,:))/norm(w_true(19:27,:),2));
        fprintf('w4 = %f, E(C_K) = %f, E(Gamma_k) = %f\n',w4(8),norm(Cmon(supp(4),tt_start_ind:tt_end_ind)-Ctrue(supp(4),tt_start_ind:tt_end_ind))...
            /norm(Ctrue(supp(4),tt_start_ind:tt_end_ind),2),norm(w4-w_true(28:36,:))/norm(w_true(28:36,:),2));
    else
        disp('Wrong equation name');
    end
else
    if isequal('cos(t)', ft_name)
        MaxIt_2 = 10;
        gamma1 = 0.0000001;
        lambda1 = 10;
        w1 = STRidge(Psi,y1,gamma1,lambda1,MaxIt_2);
        gamma2 = 0.0000001;
        lambda2 = 2;
        w2 = STRidge(Psi,y2,gamma2,lambda2,MaxIt_2);
        gamma3 = 0.000001;
        lambda3 = 0.003;
        w3 = STRidge(Psi,y3,gamma3,lambda3,MaxIt_2);
        gamma4 = 0.000001;
        lambda4 = 0.001;
        w4 = STRidge(Psi,y4,gamma4,lambda4,MaxIt_2);

        w = zeros(16,1);
        w(1:4,:) = w1;
        w(5:8,:) = w2;
        w(9:12,:) = w3;
        w(13:16,:) = w4;
        w_true = zeros(16,1);
        w_true(1,:) = alpha;
        w_true(5,:) = -1;
        w_true(9,:) = -2;
        w_true(11,:) = 0.5;
        w_true(13,:) = 2;
        w_true(15,:) = -0.5;

        %RMS Error
        fprintf('w1 = %f, E(C_K) = %f, E(Gamma_k) = %f\n',w1(1),norm(Cmon(supp(1),tt_start_ind:tt_end_ind)-Ctrue(supp(1),tt_start_ind:tt_end_ind))...
            /norm(Ctrue(supp(1),tt_start_ind:tt_end_ind),2),norm(w1-w_true(1:4,:))/norm(w_true(1:4,:),2));
        fprintf('w2 = %f, E(C_K) = %f, E(Gamma_k) = %f\n',w2(1),norm(Cmon(supp(2),tt_start_ind:tt_end_ind)-Ctrue(supp(2),tt_start_ind:tt_end_ind))...
            /norm(Ctrue(supp(2),tt_start_ind:tt_end_ind),2),norm(w2-w_true(5:8,:))/norm(w_true(5:8,:),2));
        fprintf('w3 = %f,%f, E(C_K) = %f, E(Gamma_k) = %f\n',w3(1),w3(3),norm(Cmon(supp(3),tt_start_ind:tt_end_ind)-Ctrue(supp(3),tt_start_ind:tt_end_ind))...
            /norm(Ctrue(supp(3),tt_start_ind:tt_end_ind),2),norm(w3-w_true(9:12,:))/norm(w_true(9:12,:),2));
        fprintf('w4 = %f, %f, E(C_K) = %f, E(Gamma_k) = %f\n',w4(1),w4(3),norm(Cmon(supp(4),tt_start_ind:tt_end_ind)-Ctrue(supp(4),tt_start_ind:tt_end_ind))...
            /norm(Ctrue(supp(4),tt_start_ind:tt_end_ind),2),norm(w4-w_true(13:16,:))/norm(w_true(13:16,:),2));

    elseif isequal('exp(t)', ft_name)
        MaxIt_2 = 10;
        gamma1 = 0.0000001;
        lambda1 = 10;
        w1 = STRidge(Psi,y1,gamma1,lambda1,MaxIt_2);
        gamma2 = 0.0000001;
        lambda2 = 2;
        w2 = STRidge(Psi,y2,gamma2,lambda2,MaxIt_2);
        gamma3 = 0.000001;
        lambda3 = 0.008;
        w3 = STRidge(Psi,y3,gamma3,lambda3,MaxIt_2);
        gamma4 = 0.000001;
        lambda4 = 0.008;
        w4 = STRidge(Psi,y4,gamma4,lambda4,MaxIt_2);

        w = zeros(16,1);
        w(1:4,:) = w1;
        w(5:8,:) = w2;
        w(9:12,:) = w3;
        w(13:16,:) = w4;
        w_true = zeros(16,1);
        w_true(1,:) = alpha;
        w_true(5,:) = -1;
        w_true(9,:) = 3;
        w_true(10,:) = 3;
        w_true(11,:) = 1.5;
        w_true(13,:) = -3;
        w_true(14,:) = -3;
        w_true(15,:) = -1.5;

        %RMS Error
        fprintf('w1 = %f, E(C_K) = %f, E(Gamma_k) = %f\n',w1(1),norm(Cmon(supp(1),tt_start_ind:tt_end_ind)-Ctrue(supp(1),tt_start_ind:tt_end_ind))...
            /norm(Ctrue(supp(1),tt_start_ind:tt_end_ind),2),norm(w1-w_true(1:4,:))/norm(w_true(1:4,:),2));
        fprintf('w2 = %f, E(C_K) = %f, E(Gamma_k) = %f\n',w2(1),norm(Cmon(supp(2),tt_start_ind:tt_end_ind)-Ctrue(supp(2),tt_start_ind:tt_end_ind))...
            /norm(Ctrue(supp(2),tt_start_ind:tt_end_ind),2),norm(w2-w_true(5:8,:))/norm(w_true(5:8,:),2));
        fprintf('w3 = %f,%f,%f,%f, E(C_K) = %f, E(Gamma_k) = %f\n',w3(1),w3(2),w3(3),w3(4),norm(Cmon(supp(3),tt_start_ind:tt_end_ind)-Ctrue(supp(3),tt_start_ind:tt_end_ind))...
            /norm(Ctrue(supp(3),tt_start_ind:tt_end_ind),2),norm(w3-w_true(9:12,:))/norm(w_true(9:12,:),2));
        fprintf('w4 = %f,%f,%f,%f, E(C_K) = %f, E(Gamma_k) = %f\n',w4(1),w4(2),w4(3),w3(4),norm(Cmon(supp(4),tt_start_ind:tt_end_ind)-Ctrue(supp(4),tt_start_ind:tt_end_ind))...
            /norm(Ctrue(supp(4),tt_start_ind:tt_end_ind),2),norm(w4-w_true(13:16,:))/norm(w_true(13:16,:),2));
    else
        disp('Wrong equation name');
    end
end