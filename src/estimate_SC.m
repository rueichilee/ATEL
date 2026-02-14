function [Y_hat, w, v] = estimate_SC(y, X, T0)
% ESTIMATE_SC Synthetic Control Method (Abadie, Diamond, Hainmueller 2010)
% This code is adapted from the  Synth package and tailored for our dataset.
% Description:
%   Estimates the synthetic control weights (W) and predictor weights (V).
%   The estimator minimizes the distance between the treated unit and the
%   convex combination of control units in the pre-treatment period.
%
% Inputs:
%   y  : (N+1 x T) Outcome matrix. Row 1 is Treated, Rows 2:N+1 are Controls.
%   X  : (N x T x P) Covariates matrix.
%   T0 : Integer. Length of the pre-treatment period.
%
% Outputs:
%   Y_hat : (1 x T) Vector of the estimated synthetic control outcome.
%   w     : ((N-1) x 1) Vector of weights for control units.
%   v     : (K x 1) Vector of importance weights for predictors
%% to estimate synthetic control model by Abadie, Diamond, and Hainmueller(2010)
%% This is the code original from Synth Package that is tailored for our dataset.

%% synthetic control
X = [squeeze(mean(X,2)),y(:,1:T0) ];
X = X';
X0 = X(:,2:end);
X1 = X(:,1);
% Normalization (probably could be done more elegantly)
bigdata = [X0,X1];
divisor = std(bigdata');
scamatrix = (bigdata' * diag(( 1./(divisor) * eye(size(bigdata,1))) ))';
X0sca = scamatrix([1:size(X0,1)],[1:size(X0,2)]);
X1sca = scamatrix(1:size(X1,1),[size(scamatrix,2)]);
X0 = X0sca;
X1 = X1sca;
clear divisor X0sca X1sca scamatrix bigdata;


% Y0 : control outcomes
Y0 =y(2:end,:)';
% Y1 : treated outcomes
Y1 =y(1,:)';

% Now pick Z matrices, i.e. the pretreatment period
% over which the loss function should be minmized
 
Z0 = Y0(1:T0,:);
Z1 = Y1(1:T0,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now we implement Optimization

% Check and maybe adjust optimization settings if necessary
options = optimset('fmincon');
options = optimset(options, 'Display', 'off', 'Algorithm', 'sqp');
% Get Starting Values 
s = std([X1 X0]')';
s2 = s; s2(1)=[];
s1 = s(1);
v20 =((s1./s2).^2);

lb_v = zeros(size(v20)); 

% Run Optimization (suppressed output)
[v2,fminv,exitflag] = fmincon('loss_function_sc',v20,[],[],[],[],...
   lb_v,[],[],options,X1,X0,Z1,Z0);
% display(sprintf('%15.4f',fminv));
v = [1;v2];
% V-weights

% Now recover W-weights
D = diag(v);
H = X0'*D*X0;
H = (H+H')/2;
f = - X1'*D*X0;
options = optimset('quadprog');
options = optimset(options, 'Display', 'off');
[w,fval,e]=quadprog(H,f,[],[],ones(1,size(X0,2)),1,zeros(size(X0,2),1),ones(size(X0,2),1),[],options);
w= abs(w); 
Y_hat = w' * Y0';