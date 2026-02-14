%% Process Raw Data for ATEL
% Reads 'raw_data.csv' and converts it into the .mat format 
% required by the ATEL package.

clear; clc;

% Load Raw Data
csv_file = 'raw_data.csv';
if ~exist(csv_file, 'file')
    error('File "%s" not found in current directory.', csv_file);
end

data = readtable(csv_file);

ratvio = data.ratvio;
fopsstat = data.fipsstat;

Arizona = 4; % Arizona (Treated Unit)
treated_id = fopsstat == Arizona;
control_id = fopsstat ~= Arizona;

% N: Total units (15), T: Total time periods (30)
% T0: Pre-treatment periods (17), T1: Post-treatment periods (13)
N=15; T=30; T0=17; T1=13;

% Identify unique control unit IDs (excluding Arizona)
control = unique(fopsstat,'stable');
control = control(2:end); 

y = zeros(N,T);
povrate = zeros(N,T);
policerate = zeros(N,T);

% Assign Treated Unit Data to Row 1
y(1,:) = ratvio(treated_id)';
povrate(1,:) = data.povrate(treated_id)';
policerate(1,:) = data.lnlpolicerate(treated_id)';

% Loop through Control Units to fill Rows 2 to N
for i = 1:length(control)
    % Find indices for the current control state
    idx = fopsstat == control(i);
    
    y(i+1,:) = ratvio(idx)';
    povrate(i+1,:) = data.povrate(idx)';
    policerate(i+1,:) = data.lnlpolicerate(idx)';
end

% Save Processed Data to data.mat file
save("data.mat","y","povrate","policerate");