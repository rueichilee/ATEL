function Sigma_t = cal_pointvar(F, beta, u_hat, h, T1)
% CAL_POINTVAR Calculates pointwise variance for the counterfactual path
%
% This function estimates the variance of the counterfactual
% outcome at each post-treatment time t. This is used to construct
% pointwise confidence intervals 
%
% Inputs:
%   F     : (T x J) Estimated factors.
%   beta  : (J x T1) estimated loadings for post-treatment periods.
%   u_hat : (1 x T0) Residual for the pre-treatment period.
%   h     : Scalar bandwidth.
%   T1    : Integer. Number of post-treatment periods.
%
% Output:
%   Sigma_t : (T1 x 1) Vector of variances for t = T0+1 to T.
    [T,~] = size(F);
    T0 = T-T1;

    F0 = F(1:T0,:);
    F1 = F(T0+1:end,:);
    T1h = floor(T1 * h);
    T0h = floor(T0 * h);

    KT0 = ((1:T0)-T0 ) /T0h;
    kT0 = 0.75 * (1 - KT0.^2).*(KT0 >= -1).*2 ;
    KT1 = ((T0+1:T)-T0 ) /T1h;
    kT1 = 0.75 * (1 - KT1.^2).*(KT1 <= 1).*2;

    beta_demean = beta - sum(beta .* kT1,2)/T1h;
    beta_demean = beta_demean .* sqrt(kT1);

    Sigma_t1 =  diag(F1 * beta_demean * beta_demean'* F1') ./T1h;

    xi_inv = pinv(F0'*diag(kT0)*F0);
    Sigma_t2 =  diag(F1 * xi_inv * (F0'*diag(kT0)*diag(diag(u_hat' * u_hat )) *diag(kT0) * F0) *xi_inv * F1')./T0h;

    Sigma_t = Sigma_t1 + Sigma_t2;
end
