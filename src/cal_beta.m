function beta = cal_beta(y, F0, h)
% Estimates post-treatment factor loadings 
%
% Implements Step 3 of the estimation procedure to estimate post-treatment factor loadings
%
% Inputs:
%   y  : (1 x T) Outcomes for the treated unit.
%   F0 : (T0 x J) Estimated factors for the pre-treatment period.
%   T  : Scalar. Total number of time periods.
%   h  : Scalar. Bandwidth.
%
% Output:
%   beta : (J x T1) estimated loadings for post-treatment periods.
    T = length(y);
    [T0,J] = size(F0);
    T1 = T - T0;
    T0h = floor(T0*h);
    
    %calculate kernel
    KT = ((1:T0)-T0 ) /T0h;
    kT = 0.75 * (1 - KT.^2).*(KT >= -1).*2 ;
    
    % Local Linear Estimator
    X0 = [F0,F0.*((1:T0)./T0-1)'];
    beta_T0 = pinv(X0'*diag(kT)*X0) * (X0'*diag(kT)*y(1,1:T0)'  );
    beta = beta_T0(1:J) + beta_T0(J+1:end) .* ((T0-(T0+1:T))/T1);
end

