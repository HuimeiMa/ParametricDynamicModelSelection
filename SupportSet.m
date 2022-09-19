function supp = SupportSet(L, indc, eq_name)

% ============================================================
% Description:
%   The function finds the support set of equation coefficients, 
%   such as Burger's equation and Lorenz 96 equation.
%
% Inputs:
%   L =  legend of D
%   indc = column index (restriction of the data)
%   eq_name = The name of the equation
% Outputs:
%   supp = the support set
%
% Authors: Huimei Ma, Xiaofan Lu, Linan Zhang.
% ============================================================

if isequal('lorenz96', eq_name)
    % Find the non-zeros terms of Lorenz 96 equation coefficients.
    str_u1 = strcat('u', num2str(indc(1))); % u_{1}
    str_u2 = strcat('u', num2str(indc(2))); % u_{2}
    str_un = strcat('u', num2str(indc(end))); % u_{n}
    str_unm = strcat('u', num2str(indc(end-1))); % u_{n-1}
    str_u2un = strcat(str_u2,str_un); % u_{2} * u_{n}
    str_unmun = strcat(str_unm,str_un); % u_{n-1} * u_{n}

    % Find the support setof Lorenz 96 equation coefficients.
    supp = [1,                          find(ismember(L,str_u1)), ...
        find(ismember(L,str_u2un)), find(ismember(L,str_unmun))];

else
    disp('Wrong equation name'); 
end
end
