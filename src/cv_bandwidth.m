function h_opt = cv_bandwidth(y, F)
% CV_BANDWIDTH Selects optimal bandwidth via Leave-One-Out Cross Validation
%
% This function implements the bandwidth selection procedure described in
% Eq (20). It validates the boundary estimator's ability to predict 
% pre-treatment data points.
%
% Inputs:
%   y     - (1 x T) vector of outcomes for the treated unit
%   F     - (T x J) matrix of estimated factors
%
% Output:
%   h_opt - Scalar optimal bandwidth
    hs = 0.3:0.05:0.95;
    T = length(y);
    X0 = [F,F.*((1:T)./T-1)'];    
    ssrs = zeros(size(hs,2),1);
    for ii = 1:size(hs,2)
        ssr = 0;
        h = hs(ii);
        Th = floor(T*h);
        for t = 1:T
            KT = ((1:T)-T ) /Th;
            kT = 0.75 * (1 - KT.^2).*(KT >= -1).*2 ;
            kT(t) = 0;
            beta =pinv(X0'*diag(kT)*X0) * (X0'*diag(kT)*y'  );
            m_hat = X0 * beta;
            ssr = ssr + (y(1,t) -  m_hat(t))^2;
        end
        ssrs(ii) = ssr/T;
    end
    [~,idx] = min(ssrs);
    h_opt = hs(idx);
end