
% By Jushan Bai, March 2005, New York University
% Purpose: Estimating panel data models with interactive effects

% simulation program
%  Y dependent variable , T by N 
% X independent variable  X  is T x N x p   (three dimension matrix)
% r is the true number of true factors
% F is T by r matrix of true factors
% Lambda N by r is the true loading matrix


 
% randn('state',9)

clear;
maxiter=1;  % maximum number of iterations  (set to 1000) 
Case=7;  % number of T and N combinations
BETA=zeros(maxiter, 10,Case); % contains the ``within", infeasible, the Interative effects, and the naive estimators
SIGMA=zeros(maxiter,Case);
NNN=zeros(maxiter,Case);   % the number of iterations to achieve convergence 
grandmean=zeros(maxiter,1);
p=5;    % number of regressors 
TT=[10,20,50,100,100,100 100];  % sample size T
NN=[100,100,100,100, 10,20,50];  % number of cross-sections N 
BETA2=zeros(maxiter, 10, Case);

for kk=1:Case;
    kk
    T=TT(kk);
    N=NN(kk);
    oneT=ones(T,1);
    oneN=ones(N,1);
    M_T=eye(T)-oneT*oneT'/T;   % projection matrix
    M_N=eye(N)-oneN*oneN'/N;   % projection matrix 
    
   for  niter=1:maxiter
    r=2;
    p=5;
% Data generating process
lambda=randn(N,r);   % loadings
% lambda=[randn(N,1), ones(N,1)];         % additive effect
F=randn(T,r);      % factors
% F=[ ones(T,1), randn(T,1)];             % additive effect
%e=diag(hetero(1:T))* randn(T,N);   % heteroskedastic disturbances 
e= 2*randn(T,N);   % disturbances 
XX=zeros(T,N,p);    % regressor matrix, must be T by N  by  p 

X1=randn(T,N)+ F*lambda'+ones(T,r)*lambda'+F*ones(r,N)+1;
X2=randn(T,N)+ F*lambda'+ones(T,r)*lambda'+F*ones(r,N)+1;
X3=zeros(T,N);   %initialization
X4=X3;           % initialization
X5=X3;           % initialization
X3(1:T,1:N)=1;   % grand mean
aa=sum(lambda')+randn(1,N);
bb=sum(F,2)+randn(T,1);
%aa=lambda(:,1).^3'-randn(1,N).^2*0.1;
%bb=F(:,2).^3-randn(T,1).^2*0.1;

for t=1:T;
    X4(t,:)=aa;
end
for i=1:N;
    X5(:,i)=bb;
end

truebeta=[1,3,5,2,4];
XX(:,:,1)=X1;
XX(:,:,2)=X2;
XX(:,:,3)=X3;
XX(:,:,4)=X4;
XX(:,:,5)=X5;
X=XX(:,:,1:p);
Y=F*lambda'+ e;
for k=1:p;
    Y=Y+X(:,:,k)*truebeta(k);
end         % Data generating is done

%res=Y-X1-X2*3-X4*2-X5*4;
%res1=res-F*lambda'-e;
%mean(mean(res1))
%grand=mean(mean(res));

%meanY=mean(mean(Y));
%meanX2=mean(mean(X2));
%Y=Y-meanY;
%X2=X2-meanX2;;

MFF=zeros(T,T,3);    % to contain three projection matrices
MLL=zeros(N,N,3);
Beta0=zeros(p,3);    % to contain three different estimates

MF=eye(T)-F*inv(F'*F)*F';   %projection matrix for the infeasible estimator

% Infeasible estimators:
% betaL and betaF are different since they are estimating different 
% number of parameters.
% ML=eye(N)-lambda*inv(lambda'*lambda)*lambda';  % projection matrix 
% [betaL]=panelbeta(X1',X2',Y',ML)

  [betaF]=Mul_panelbeta(X,Y,MF) ;   % the infeasible estimator   
      
  % Within transformation:   (not used in this version)
    X1FL=M_T*X1*M_N;
    X2FL=M_T*X2*M_N;
    YFL=M_T*Y*M_N;
  
 %  [betaWG]=panelbeta(X1FL,X2FL,YFL,eye(T));   % the within estimator

    [beta0]=Mul_panelbeta(X,Y,eye(T));  % naive estimator (also used as a starting value)
    BETA2(niter, 1:p, kk)=beta0';    % ols estimator
    U=Y;
    for k=1:p;
        U=U-X(:,:,k)*beta0(k);
    end
    
    [F1,L1,VNT]=panelFactorNew(U,r); 
    
    
    
   %[beta0]=panelbeta(MF*X1,MF*X2,MF*Y,eye(T))   % removing individual fixed effects
   %[beta0]=panelbeta(X1*ML,X2*ML,Y*ML,eye(T))   % removing time effects
  % the above can be used as an initial value for beta
  %  [beta0]=panelbeta(X1FL,X2FL,YFL,eye(T));   % the within estimator (also used as a staring value)
  %  U=Y-X1*beta0(1)-X2*beta0(2);
  %  [F1,L1,VNT]=panelFactor(U,r);
  %  MFF(:,:,3)=eye(T)-F1*F1'/T;
   % Beta0(:,3)=beta0;
  
     
    beta=zeros(p,3);    % to contain the interative effect estimators, with different staring methods (3 methods)
    nnn=zeros(1,3);   % contain the number of iterations for the three methods to achieve convergence
    sigma2=zeros(1,3); % the estimated residual variance for each starting method, also the optimal value of the objective function
  
    % compute (X'X)^{-1}
    
     [XXinv]= Mul_XXinv(X);   % outside the beta iteration loop
          [beta(:,1), F1,L1, VNT, e1, nnn(1)]=Mul_betaIterNew(X,XXinv, Y, F1,L1, r, 0.0001);
          sigma2(1)=trace(e1*e1')/(N*T-r*(N+T)+r^2-2);
    
       BETA(niter,1:p,kk)=beta(:,1)';    % the interactive effects estimator
       NNN(niter,kk)=nnn(1);
       SIGMA(niter,kk)=sigma2(1);
    %  BETA(niter, 6:7,kk)=betaWG';      % within  estimator
       BETA(niter, (p+1):(2*p),kk)=betaF';      % infeasible
     
   
     
end  % niter loop
end  % kk loop

save('Mul_BETA', 'BETA');
save('Mul_NNN', 'NNN');
save('Mul_SIGMA', 'SIGMA');

 