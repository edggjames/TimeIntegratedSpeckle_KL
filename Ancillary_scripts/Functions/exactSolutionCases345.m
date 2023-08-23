function [pW_all,W_bar] = exactSolutionCases345(case_number,T,N,t_c,thresh,res,I_vec,alpha)
% Function to calculate exact solution PDFs from Goodman, for test cases 3,
% 4 and 5

if ~exist('alpha','var')
    % if alpha is not parsed set to default value of 1
    disp('Setting alpha to default value of 1 in Mone')
    alpha = 1;
else

num_Ts = length(T);

% loop through T vals
pW_all = zeros(num_Ts,res);
W_bar  = zeros(num_Ts,1);

for i = 1:num_Ts
    
    T_temp = T(i);
    t2 = linspace(0,T_temp,N);
    t1 = t2;
    [n,k] = meshgrid(t2,t1);
    diff = abs(k-n);
    if case_number == 3
        K = alpha*exp(-diff/t_c) + (1-alpha);
    elseif case_number == 4
        K = alpha*exp(-sqrt(diff/(t_c))) + (1-alpha);
    elseif case_number == 5
        K = alpha*exp(-(diff/(t_c)).^2)  + (1-alpha);
    end
    e = eigs(K,N)/N;
    
    %% truncate - how many eigen values to keep
    e(e<thresh) = [];
    num_lambda = length(e);
    disp(num_lambda)

    %% take sum over lambda_n
    
    pW = zeros(1,res);
    d_n_vec = zeros(num_lambda,1);
    
    for n = 1:num_lambda
        lambda_n = e(n);
        % take product over d_n
        d_n = 1;
        for m = 1:num_lambda
            lambda_m = e(m);
            if n ~= m
                multiplier = 1/(1-(lambda_m/lambda_n));
                d_n = d_n * multiplier;
            end
        end
        d_n_vec(n) = d_n;
        pW_temp = d_n/lambda_n*(exp(-I_vec/lambda_n));
        % this is the key step - each of these distributions is an image
        % of mean lambda_n multiplied by d_n
        pW = pW + pW_temp;
    end
    
    % check sum of d_n_vec
    d_n_vec_sum = sum(d_n_vec);
    disp(d_n_vec_sum)
    
    % check area under curve using the trapezoid rule
    h = (I_vec(end) - I_vec(1))/res;
    area = (h/2)*(2*sum(pW) - pW(1) - pW(end));
    disp(area)
    
    % assign
    pW_all(i,:) = pW;
    
    % calculate mean
    W_bar(i) = sum(d_n_vec.*e);
    
end
end

