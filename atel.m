function results = atel(Y, X, T0, varargin)
% ATEL Average Treatment Effect Localization Estimator
%
%   RESULTS = ATEL(Y, X, T0) estimates the Short-Term Average Treatment Effect
%   using Diversified Projection (DP) and Local Linear Smoothing.
%
%   INPUTS:
%       Y  : (N+1) x T matrix of outcomes.
%            - Row 1: Treated unit.
%            - Rows 2 to N+1: Control units.
%       X  : (N+1) x T x P array of observed covariates.
%       T0 : Scalar integer. The last period of the pre-treatment window.
%
%   OPTIONAL PARAMETERS (Name-Value pairs):
%       'Rank'      : Integer (J). Number of factors/sieve terms. Default: 2.
%       'Basis'     : String. Basis for weights W. 
%                     Options: 'bspline' (default), 'polynomial', 'trigonometric'.
%       'Bandwidth' : Scalar (h). Bandwidth for kernel smoothing. 
%                     Default: 'auto' (calculated via Cross-Validation).
%       'Alpha'     : Scalar. Significance level (e.g., 0.05 for 95% CI).
%
%   OUTPUTS (structure):
%       results.atel      : Scalar point estimate of ATEL.
%       results.se        : Standard error.
%       results.ci        : [Lower, Upper] Confidence Interval.
%       results.alpha_t   : Vector of treatment effects over time.
%       results.Y1_hat    : Vector of counterfactual outcomes.
%       results.se_point  : Vector of standard error for counterfactual outcomes
%       results.optimal_h : Scalar of optimal bandwidth
%
%   REFERENCES:
%       Lee, R.-C. (2026). "Average Treatment Effect Localization: Projection
%       Methods in Synthetic Control".

    %% 1. Input Parsing & Setup
    p = inputParser;
    addRequired(p, 'Y', @isnumeric);
    addRequired(p, 'X', @isnumeric);
    addRequired(p, 'T0', @isscalar);
    addParameter(p, 'Rank', 2, @(x) isscalar(x) && x > 0);
    addParameter(p, 'Basis', 'bspline', @(x) any(validatestring(x, {'bspline', 'polynomial', 'trigonometric'})));
    addParameter(p, 'Bandwidth', [], @(x) isempty(x) || (isscalar(x) && x > 0));
    addParameter(p, 'Alpha', 0.05, @(x) isscalar(x) && x > 0 && x < 1);
    
    parse(p, Y, X, T0, varargin{:});
    
    J = p.Results.Rank;
    basis_type = p.Results.Basis;
    sig_level = p.Results.Alpha;
    
    [N_total, T] = size(Y);
    N = N_total - 1;            % Number of control units
    T1 = T - T0;                % Number of post-treatment periods
    
    % Separation of Treated (Unit 1) and Control (Units 2:N+1)
    Y1 = Y(1, :);               
    Y_N = Y(2:end, :);          
    
    %% 2. Step 1: Construct Diversified Weights (W)
    W = construct_weights(X, J, basis_type);

    %% 3. Step 2: Estimate Factors (F_t) via Diversified Projection
    F = zeros(T,J);
    for j=1:J
        F(:,j) = mean(Y_N .* W(2:end,1+(j-1)*T:j*T),1)';
    end

    F0 = F(1:T0,:);
    F1 = F(T0+1:end,:);
    
    %% 4. Step 3: Estimate Factor Loadings (beta_1t)
    
    % Bandwidth Selection (Cross-Validation)
    if isempty(p.Results.Bandwidth)
        h = cv_bandwidth(Y1(1:T0), F0);
    else
        h = p.Results.Bandwidth;
    end
    T0h = floor(T0*h);
    T1h = floor(T1*h);

    beta = cal_beta(Y1, F0, h);
    %% 5. Step 4: Construct Counterfactuals 
    Y1_hat_post = sum(beta .* F1',1);

    %% 6. Step 5: Calculate ATEL 
    Y1_post = Y1(T0+1:end);
    alpha_t = Y1_post - Y1_hat_post; 
    
    % Calculate kernel weights
    KT1 = ((T0+1:T)-T0 ) /T1h;
    kT1 = 0.75 * (1 - KT1.^2).*(KT1 <= 1).*2;

    atel_val = sum(alpha_t.*kT1) / T1h  ;
    
    %% 7. Inference (Asymptotic Normality)

    % Calculate standard error
    [Sigma, Y1_hat_pre] = cal_var(Y1,F0,F1,h);
    se_atel = sqrt(Sigma);

    % construct counterfactul outcomes
    Y1_hat = [Y1_hat_pre, Y1_hat_post];

    % confidence interval
    crit_val = norminv(1 - sig_level/2);
    ci_lower = atel_val - crit_val * se_atel;
    ci_upper = atel_val + crit_val * se_atel;

    % Calculate P-value
    t_stat = atel_val / se_atel;
    p_val = 2 * tcdf(-abs(t_stat), T1h-1);

    % Calculate pointwise standard error for counterfactual outcomes
    Sigma_point = cal_pointvar(Y_N,F,h,T1);
    se_point = sqrt(Sigma_point);
    %% 8. Pack Results
    results.atel = atel_val;
    results.se = se_atel;
    results.ci = [ci_lower, ci_upper];
    results.alpha_t = alpha_t;
    results.Y1_hat = Y1_hat;
    results.optimal_h = h;
    results.p_val = p_val;
    results.se_point=se_point;
    results.W = W;
    results.F = F;
    results.beta = beta;
end