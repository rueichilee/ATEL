function [k_IC1, k_IC2, IC1, IC2] = select_rank(y, X, basis_type, kmax)
% IC determined the number of rank:  Su & Wang (2017) On time-varying
%  factor models: Estimation and testing.
%
%   INPUT:
%       Y  : (N+1) x T matrix of outcomes.
%            - Row 1: Treated unit.
%            - Rows 2 to N+1: Control units.
%       X  : (N+1) x T x P array of observed covariates.
%       Js : desired basis dimension (2, 3, 4, or 5)
%       basis_type : String: 'bspline', 'polynomial', or 'trigonometric'
%       kmax : (optional) maximum number of candidate factors.
%              Default: min(8, m-1) where m = min(N, T).
%
%   OUTPUT:
%       k_IC1  : estimated number of factors using the IC1 criterion
%       k_IC2  : estimated number of factors using the IC2 criterion
%       IC1    : kmax-by-1 vector of IC1
%       IC2    : kmax-by-1 vector of IC2
%   Reference: Ahn and Horenstein (2013), "Eigenvalue Ratio Test for
%   the Number of Factors", Econometrica.

    % Dimensions
    [N_total, T] = size(y);
    N = N_total - 1;            % Number of control units

    m = min(N, T);

    % Default kmax
    if nargin < 2 || isempty(kmax)
        kmax = min(8, m - 1);  % you can adjust this rule if desired
    else
        kmax = min(kmax, m - 1);
    end

    %each sereis is demeand and standardized 
    y  = y - mean(y,2);             % demean over T
    ys  = std(y, 0, 2);
    ys(ys< 1e-12) = 1;              % avoid divide-by-zero
    y  = y ./ ys;

    X = X - mean(X, 2);  % demean over T 
    xs = std(X, 0, 2);   
    xs(xs < 1e-12) = 1;             % avoid divide-by-zero
    X  = X ./ xs;                   

    IC1 = zeros(kmax-1, 1);
    IC2 = zeros(kmax-1, 1);

    for J = 2:kmax
        W = construct_weights(X, J, basis_type);
        F = zeros(T,J);
        for j=1:J
            F(:,j) = mean(y(2:end,:) .* W(2:end,1+(j-1)*T:j*T),1)';
        end

        %select bandwidth
        hs = 0.3:0.05:0.95;
        X0 = [F,F.*((1:T)./T-1)'];    
        ssrs = zeros(size(hs,2),1);
        for ii = 1:size(hs,2)
            ssr = 0;
            h = hs(ii);
            for t = 1:T
                Th = floor(T*h);
                KT = ((1:T)-T ) /Th;
                kT = 0.75 * (1 - KT.^2).*(KT >= -1).*2 ;
                kT(t) = 0;
                beta =pinv(X0'*diag(kT)*X0) * (X0'*diag(kT)*y(2:end,:)'  );
                fitted = X0 * beta;
                ssr = ssr + sum((y(2:end,t) -  fitted(t,:)').^2);
            end
            ssrs(ii) = ssr/(N*T);
        end
        [~,idx] = min(ssrs);
        h = hs(idx);
        
        % cross validation
        Th = floor(T*h);
        V=0;
        for t = 1:T
            K = ((1:T) - t ) / Th;
            if t < Th
                tmp = t/Th;
                tmp = (tmp + (tmp^3)/3 + 2/3) * 3/4; % denomiator adjustment
                k = 0.75 * (1 - K.^2) .* (abs(K)<=1) ./ tmp;
            elseif t > T - Th
                tmp = (1-t/T)/h;
                tmp = (tmp - (tmp^3)/3 + 2/3) * 3/4; % denomiator adjustment
                k = 0.75 * (1 - K.^2) .* (abs(K)<=1) ./ tmp;
            else
                k = 0.75 * (1 - K.^2) .* (abs(K)<=1);
            end
            X00 = [F,F.*(((1:T) - t) / T)'];
            beta = pinv(X00'*diag(k)*X00) * (X00'*diag(k)*y(2:end,1:T)'  );
            F_tilde = pinv(beta*beta') * beta*y(2:end,t);
            tmp = y(2:end,t) - beta' * F_tilde;
            V = V + sum(tmp.^2);
        end
        V = V/(N*T);
        
        % The penalty factor common to both (approximates complexity cost)
        % Note: (N+T)/(N*T) is equivalent to (1/T + 1/N)
        sparsity_factor = (N + Th) / (N * Th);
        
        % --- IC_p1 ---
        % Penalty uses log of the average sample size
        penalty1 = J * sparsity_factor * log((N * Th) / (N + Th));
        IC1(J-1) = log(V) + penalty1;
        
        % --- IC_p2 ---
        % Penalty uses log of the minimum dimension (more conservative)
        penalty2 = J * sparsity_factor * log(min(N, Th));
        IC2(J-1) = log(V) + penalty2;
    end
    [~, k_IC1] = max(IC1);
    [~, k_IC2] = max(IC2);
    % J start with 2
    k_IC1 = k_IC1 + 1; 
    k_IC2 = k_IC2 + 1;
end