function [Sigma, pre_Y_hat] = cal_var(y,F0,F1,h)
% CAL_VAR Calculates the asymptotic variance of the ATEL estimator
%
% This function estimates the variance components defined in Theorem 1.
% It consists of two parts:
%   Sigma1 : Variance arising from the estimation error of factor loadings (beta).
%   Sigma2 : Variance arising from the idiosyncratic errors (u_t).
%
% Inputs:
%   y  : (1 x T) Outcome vector for the treated unit.
%   F0 : (T0 x J) Estimated factors for the pre-treatment period.
%   F1 : (T1 x J) Estimated factors for the post-treatment period.
%   h  : Scalar bandwidth.
%
% Output:
%   Sigma : scalar variance.
%   Pre_Y_hat : (1 x T0) Estimated Y hat for the pre-treatment period.
    [T0,J] = size(F0);
    T = length(y);
    T1 = T-T0;
    T0h = floor(T0*h);
    T1h = floor(T1*h);
    
    %calculate residual and hat Y1 in the pre-treatment period
    pre_Y_hat = zeros(1,T0);
    u_hat = zeros(1,T0);
    for t = 1:T0
        K = ((1:T0) - t ) / T0h;
        if t < T0h
            tmp = t/T0h;
            tmp = (tmp + (tmp^3)/3 + 2/3) * 3/4; % denomiator adjustment
            k = 0.75 * (1 - K.^2) .* (abs(K)<=1) ./ tmp;
        elseif t > T0 - T0h
            tmp = (1-t/T0)/h;
            tmp = (tmp - (tmp^3)/3 + 2/3) * 3/4; % denomiator adjustment
            k = 0.75 * (1 - K.^2) .* (abs(K)<=1) ./ tmp;
        else
            k = 0.75 * (1 - K.^2) .* (abs(K)<=1);
        end
        X00 = [F0,F0.*(((1:T0) - t) / T0)'];
        beta = pinv(X00'*diag(k)*X00) * (X00'*diag(k)*y(1,1:T0)'  );
        fitted = X00 * beta;
        pre_Y_hat(t) = fitted(t);
        u_hat(t) = y(1,t) - pre_Y_hat(t);
    end

    %calculate kernel
    KT0 = ((1:T0)-T0 ) /T0h;
    kT0 = 0.75 * (1 - KT0.^2).*(KT0 >= -1).*2 ;
    KT1 = ((T0+1:T)-T0 ) /T1h;
    kT1 = 0.75 * (1 - KT1.^2).*(KT1 <= 1).*2;

    %calculate Sigma1
    lambda =  sum(F1'.*kT1,2) / T1h;
    xi_inv = pinv(F0'*diag(kT0)*F0);
    Sigma1 = lambda' * xi_inv * (F0'*diag(kT0)*diag(diag(u_hat' * u_hat )) *diag(kT0) * F0) *xi_inv*lambda;

    %calculate Sigma2
    u0_hat = u_hat(1,T0-T0h+1:T0);
    KT0_ = ((T0-T0h+1:T0)-T0 ) /T0h;
    kT0_ = 0.75 * (1 - KT0_.^2).*(KT0_ >= -1).*2 ;
    Sigma2 = mean( (u0_hat.^2) .* (kT0_.^2))/T1h;

    Sigma = Sigma1 + Sigma2;
end