%% ATEL Simulation Test
% Description: Simulates data from a factor model with a KNOWN treatment 
% effect to verify that the estimator recovers the truth.

clear; clc; close all;

rng(1);
fprintf('Running ATEL Simulation Test...\n');

%% 1. Simulation Setup (DGP)
% We create data where Y = X*theta + L*F' + alpha + u
N = 50;       % Number of units
T = 60;       % Number of time periods
T0 = 40;      % Number of pre-treatment periods
T1 = T - T0;  % Number of post-treatment periods
P = 1;        % Number of covariates
J_true = 2;   % True number of factors

% Generate Latent Factors (F) and Loadings (L)
L1 = 1/sqrt(2) + randn(N,T);
L2 = 1/sqrt(2) + randn(N,T);
F = 1/sqrt(2) + randn(J_true,T);

% Generate Covariates (X) and theta
X = randn(N, T, P);
true_theta = randn(P,T);

h=0;
for p = 1:P
    h = h + squeeze(X(:,:,p)) .* true_theta(p,:) + L1.* F(1,:) + L2.* F(2,:);
end

% Construct Baseline Outcome Y(0)
u = randn(N,T);  
Y0 = h + u;

% Add Treatment Effect (Tau) to Unit 1 ONLY
alpha = 5.0; % True effect size
Y = Y0;
Y(1, T0+1:end) = Y(1, T0+1:end) + alpha;


%% 2. Run Estimator
% Note: We feed the raw Y and X. The function handles the rest.

Rank = J_true; % fix Rank=J_true for this test.
basis = 'bspline'; % 'bspline', 'polynomial', or 'trigonometric'
alpha = 0.5; % Significance level (e.g., 0.05 for 95% CI).

results = atel(y, X, T0, ...
               'Rank', Rank, ...
               'Basis', basis, ...
               'Alpha', alpha, ...
               'Bandwidth', []); % Auto-select bandwidth via CV

% Calculate kernel weights
T1h = floor(T1 * results.optimal_h);
KT1 = ((T0+1:T)-T0 ) /T1h;
kT1 = 0.75 * (1 - KT1.^2).*(KT1 <= 1).*2;
% Calculate true atel
true_atel = sum(alpha.*kT1) / T1h  ;


%% 3. Output (example)

fprintf('\n-----------------------------------------\n');
fprintf('True Effect:      %.4f\n', true_atel);    % true ATEL
fprintf('Estimated (ATEL): %.4f\n', results.atel); % scalar ATEL
fprintf('95%% CI:              [%.4f, %.4f]\n', results.ci(1), results.ci(2)); % (1-alpha)% CI
fprintf('P-value:             %.4f\n', results.p_val); % p_value under alpha significance level
fprintf('-----------------------------------------\n');

%% 4. Visualization
years = 1:T;
figure; hold on;

% Plot Actual Treated Unit
plot(years, Y(1, :), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Observed (Treated)');

% Plot Estimated Counterfactual
plot(years, results.Y1_hat(:), 'r--', 'LineWidth', 2, 'DisplayName', 'Estimated Counterfactual');

% Plot TRUE Counterfactual (Y0) - ideally they should overlap
plot(years, Y0(1, :), ':', 'LineWidth', 2, 'DisplayName', 'True Counterfactual');

xline(T0+1, '--k', 'Treatment Start', 'HandleVisibility', 'off');
legend('Location', 'best');
title('Simulation');

hold off