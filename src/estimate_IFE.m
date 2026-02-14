function results = estimate_IFE(y,X,T0,J,h)
% ESTIMATE_IFE Interactive Fixed Effects (Bai, 2009)

% Description:
%   Estimates the counterfactual using the Interactive Fixed Effects (IFE) model.
%   It uses control units to estimate common factors (F) and slope coefficients (beta),
%   then projects the treated unit's pre-treatment data onto F to get specific loadings.
%
% Inputs:
%   X  : N x T x P covariate matrix (Requires permutation for IFE)
%   y  : (N+1) x T outcome matrix (Row 1 = Treated)
%   J  : Number of factors
%   T0 : Pre-treatment period length
%   h  : Bandwidth for ATEL weighting
%
% Output (structure):
%   results.atel      : Scalar point estimate of ATEL.
%   results.se        : Standard error.
%   results.alpha_t   : Vector of treatment effects over time.
%   results.Y1_hat    : Vector of counterfactual outcomes.

    [N_total,T] = size(y);
    N = N_total - 1;
    T1 = T - T0;
    T1h = T1 * h;
    T0h = T0 * h;
    %% interactive fixed effect Bai 2009
    % Transform X from (N x T x P) -> (T x N x P) for the IFE function
    XX = permute(X, [2, 1, 3]);
    pp = size(XX,3);
    
    % Use Control Units (Rows 2:end) to estimate common structure
    [beta, F] = IFE(y(2:end,:)',XX(:,2:end,:),J);
    beta = beta(end,1:J);
    
    %% Estimate Loadings for Treated Unit (Projection)
    X_treated = squeeze(XX(:, 1, :)); % T x P
    fitted_parametric = X_treated * beta'; % T x 1
    
    % Get Pre-treatment Factors
    F0 = F(1:T0, :); 
    
    % Outcome - Parametric = Factor Component
    y_treated_pre = y(1, 1:T0)';
    resid_for_projection = y_treated_pre - fitted_parametric(1:T0);
    
    % OLS Projection to get lambda_1 (beta1 in your code)
    lambda_1 = (F0' * F0) \ (F0' * resid_for_projection);
    
    %% Construct Counterfactual & ATEL
    % Y_hat = X_1*beta + F*lambda_1
    fitted_factors = F * lambda_1; 
    Y1_hat = fitted_factors' + fitted_parametric'; % 1 x T vector
    %kernel
    KT1 = ((T0+1:T)-T0 ) /T1h;
    kT1 = 0.75 * (1 - KT1.^2).*(KT1 <= 1).*2;
    
    alpha_t = y(1, T0+1:T) - Y1_hat(1,T0+1:end);
    atel_val = sum(alpha_t.*kT1) / T1h;
    %% Inference (Variance Calculation)
    %calculate residuals
    u_hat = y(1, 1:T0) - Y1_hat(1, 1:T0);
    
    %calculate variance
    hat_F = sum(F(T0+1:end,:)'.*kT1,2) / T1h;
    hat_FF_inv = inv((F0'*F0) / T0);   
    Fu = F0.*(u_hat');
    hat_v = Fu'*Fu / T0;
    Sigma1 = hat_F'*hat_FF_inv*hat_v*hat_FF_inv*hat_F /(T0);
    
    KT0 = ((1:T0)-T0 ) /T0h;
    kT0 = 0.75 * (1 - KT0.^2).*(KT0 >= -1).*2 ;
    Sigma2 = mean( (u_hat.^2) .* (kT0.^2))/T1h;
    
    Sigma = Sigma1 + Sigma2;
    Std = sqrt(Sigma);
    
    %calculate t statistics
    t_stat = atel_val / Std;
    
    % Calculate the p-value (two-tailed test, using t-distribution)
    p_value = 2 * tcdf(-abs(t_stat), T1-1);
    %% Pack Results
    results.atel = atel_val;
    results.se = Std;
    results.alpha_t = alpha_t;
    results.Y1_hat = Y1_hat;
    results.p_val = p_value;
end