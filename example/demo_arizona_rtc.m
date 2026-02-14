%% Demo: ATEL of Right-to-Carry (RTC) Laws on Violent Crime
%
% Replicates Section 5 of Lee (2026).
% Analyzes the impact of RTC laws in Arizona (implemented in 1994).
%
% Data Specs:
%   - Time: 1977 - 2006 (30 years)
%   - Treated Unit: Arizona
%   - Control Units: 14 states (9 never treated, 5 treated after 2006) [cite: 690]
%   - Outcome: Violent crime rate per 100,000 residents
%   - Covariates: Police rate, Poverty rate [cite: 691]
%
% Methodology:
%   - Method: Diversified Projection (DP) with B-spline weights
%   - Rank (J): 2
%   - Kernel: Epanechnikov

clear; clc; close all;

% Add package to path if not already added
if exist('atel', 'file') ~= 2
    % Use genpath to add src AND src/IFE automatically
    addpath(genpath('../src')); 
    addpath('../'); 
end

%% 1. Load Data
% load 'arizona_crime.mat'. 
data_path = 'data/data.mat';
if exist(data_path, 'file')
    fprintf('Loading dataset: %s\n', data_path);
    load(data_path);
    % Expected variables in .mat file:
    %   Y          : (15 x 30) Outcome matrix
    %   povrate    : (15 x 30) Poverty rate variable
    %   policerate : (15 x 30) Police rate variable
else
    error('ATEL:DataNotFound', ...
          ['Dataset not found at: %s\n\n' ...
           'Please ensure ''arizona_crime.mat'' is placed in the ''examples/data'' folder.\n' ...
           'See README for instructions on obtaining replication data.'], data_path);
end

%% 2. Configuration
% Create X
X = cat(3, povrate, policerate);

% Determine basis function
basis  = 'bspline';

% Policy implemented in 1994, so T0 is 1993.
T0 = 1993 - (1977 - 1); 

% Specifiy rank J
J_rank = 2;       

% Rank Selection (Auto-Detect if J_rank is missing)
if isempty(J_rank)
    fprintf('J_rank not specified. Performing Rank Selection...\n');

    k_max = 8;
    
    % Call the helper function select_rank.m
    [k_IC1, k_IC2, IC1, IC2] = select_rank(y, X, basis, k_max);
    
    % Conservative choice: min of IC1 and IC2
    J_rank = min(k_IC1, k_IC2);
    
    fprintf('  - IC1 suggests: k = %d\n', k_ER);
    fprintf('  - IC2 suggests:     k = %d\n', k_GR);
    fprintf('  -> SELECTION: J = %d\n', J_rank);
else
    fprintf('Using user-specified Rank: J = %d\n', J_rank);
end


alpha  = 0.05;    % 95% Confidence Interval

%% 3. Estimation
% Calls the main package function
results = atel(y, X, T0, ...
               'Rank', J_rank, ...
               'Basis', basis, ...
               'Alpha', alpha, ...
               'Bandwidth', []); % Auto-select bandwidth via CV

%% 4. Display Results (Table 3)
fprintf('\n--- Estimation Results (Table 3 Replication) ---\n');
fprintf('Number of Factors (J): %d\n', J_rank);
fprintf('ATEL Estimate:       %.4f\n', results.atel);
fprintf('Standard Error:      %.4f\n', results.se);
fprintf('95%% CI:              [%.4f, %.4f]\n', results.ci(1), results.ci(2));
fprintf('P-value:             %.4f\n', results.p_val);

%% 5. Visualization (Figure 3)
% Plot Actual vs. Counterfactual for the Treated Unit
years = (1977:2006)';          % 30 x 1
T     = numel(years);

% Logical index for post-treatment period (1994–2006: 13 years)
idx_post = (years >= 1994);

% Build error vectors: NaN before 1994, band after 1994
band = results.se_point .* 1.96;
err = NaN(T,1); err(idx_post) = band(:);

% Treated path as column
y_treated = y(1,:)';
Y_hat = results.Y1_hat;

figure; hold on

% 1) Treated unit (no CI)
h_treated = plot(years, y_treated, 'k-', 'LineWidth', 2);
% 2) Synthetic estimates with post-treatment CIs
h2 = errorbar(years, Y_hat, err, err, 'o--','LineWidth', 1.5);

% Vertical treatment line (not in legend)
xline(1994, '--k', 'HandleVisibility', 'off');

% Axis label + sizes
ylabel('Violent Crime Rate Per 100K Residents', 'FontSize', 16);

% Legend (handles are now both valid)
legend([h_treated, h2], ...
       {'Treated Unit', '95% CI'}, ...
       'Location', 'northwest', 'FontSize', 14);

% Tick label size
set(gca, 'FontSize', 14);
ylim([300 750])
hold off;

%% 6 Estimate interactive fixed effect and synthetic control
results_IFE = estimate_IFE(y,X,T0,J_rank,results.optimal_h);
[Y_hat_SC,~,~] = estimate_SC(y, X, T0);

%% 7. Display Results (Table 4)

fprintf('\n=======================================================\n');
fprintf('             Table 4: Method Comparison                \n');
fprintf('=======================================================\n');
fprintf(' %-12s | %-10s | %-10s | %-10s \n', 'Method', 'Estimate', 'Std. Err', 'p-value');
fprintf('-------------------------------------------------------\n');

% Row 1: ATEL (Proposed)
fprintf(' %-12s | %10.4f | %10.4f | %10.4f \n', ...
        'ATEL', results.atel, results.se, results.p_val);

% Row 2: IFE (Benchmark)
fprintf(' %-12s | %10.4f | %10.4f | %10.4f \n', ...
        'IFE', results_IFE.atel, results_IFE.se, results_IFE.p_val);

%% 8. Visualization (Figure 4)
% Comparison of Counterfactuals: Actual vs. SCM vs. IFE vs. DP (ATEL)
% Treated path
y_treated = y(1,:);

figure; hold on; 

% 1) Actual Data 
plot(years, y_treated, '-', 'LineWidth', 2.5, 'DisplayName', 'Treated Unit');

% 2) Synthetic Control
plot(years, Y_hat_SC, '--', 'LineWidth', 2.5, 'DisplayName', 'SCM');

% 3) IFE 
plot(years, results_IFE.Y1_hat, ':', 'LineWidth', 2.5, 'DisplayName', 'IFE');

% 4) DP 
plot(years, results.Y1_hat, '-.', 'LineWidth', 2.5, 'DisplayName', 'DP (Proposed)');

set(gca, 'FontSize', 14);
legend('Treated Unit','SCM','IFE','DP', 'FontSize', 14)
xline(1994,'--','HandleVisibility','off');
ylabel("Violent Crime Rate Per 100K Residents", 'FontSize', 16)
hold off;




