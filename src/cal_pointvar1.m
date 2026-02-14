function Sigma_t = cal_pointvar1(y,F,h,T1)
% CAL_POINTVAR Calculates pointwise variance for the counterfactual path
%
% This function estimates the variance of the counterfactual
% outcome at each post-treatment time t. This is used to construct
% pointwise confidence intervals 
%
% Inputs:
%   y  : (1 x T) Outcome vector of the control units.
%   F  : (T x J) Estimated factors.
%   h  : Scalar bandwidth.
%   T1 : Integer. Number of post-treatment periods.
%
% Output:
%   Sigma_t : (T1 x 1) Vector of variances for t = T0+1 to T.
    [N,T] = size(y); % control outcomes
    [~,J] = size(F)
    T0 = T-T1;

    F0 = F(1:T0,:);
    F1 = F(T0+1:end,:);
    T1h = floor(T1 * h);
    T0h = floor(T0 * h);

    betat = zeros(2*J,T1);
    Sigma_t = zeros(T1,1);

    KT0 = ((1:T0)-T0 ) /T0h;
    kT0 = 0.75 * (1 - KT0.^2).*(KT0 >= -1).*2 ;
    
    % Local Linear Estimator 
    X0 = [F0,F0.*((1:T0)./T0-1)'];
    beta = pinv(X0'*diag(kT0)*X0) * (X0'*diag(kT0)*y(1,1:T0)'  );
    betat = beta(1:J) - beta(J+1:end) .* (((T0+1:T)-T0 ) /T1);

    Sigma_t = diag(F1 * betat(:,1:T1h) * betat(:,1:T1h)' * F1' ./T1h);
end
