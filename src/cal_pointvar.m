function Sigma_t = cal_pointvar(y,F,h,T1)
% CAL_POINTVAR Calculates pointwise variance for the counterfactual path
%
% This function estimates the variance of the counterfactual
% outcome at each post-treatment time t. This is used to construct
% pointwise confidence intervals 
%
% Inputs:
%   y  : (N x T) Outcome vector of the control units.
%   F  : (T x J) Estimated factors.
%   h  : Scalar bandwidth.
%   T1 : Integer. Number of post-treatment periods.
%
% Output:
%   Sigma_t : (T1 x 1) Vector of variances for t = T0+1 to T.
    [N,T] = size(y); % control outcomes
    T0 = T-T1;
    Sigma_t = zeros(T1,1);
    Th = T * h;
    for t = T0+1:T
        K = ((1:T) - t ) / Th;
        if t > T - Th % right boundary t
            tmp = t/Th;
            tmp = (tmp + (tmp^3)/3 + 2/3) * 3/4; % denomiator adjustment
            k = 0.75 * (1 - K.^2) .* (abs(K)<=1) ./ tmp;
        elseif t < Th % left boundary t
            tmp = (1-t/T)/h;
            tmp = (tmp - (tmp^3)/3 + 2/3) * 3/4; % denomiator adjustment
            k = 0.75 * (1 - K.^2) .* (abs(K)<=1) ./ tmp;
        else
            k = 0.75 * (1 - K.^2) .* (abs(K)<=1);
        end
        X00 = [F,F.*(((1:T) - t) / T)'];
        beta = pinv(X00'*diag(k)*X00) * (X00'*diag(k)*y'  );
        fitted = X00 * beta;
        fitted_demean = fitted(t,:) - mean(fitted(t,:));
        var_t = sum(fitted_demean.^2) / (N-1);
        Sigma_t(t-T0) = var_t;
    end
