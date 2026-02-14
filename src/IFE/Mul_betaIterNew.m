% OUTPUT:
% betanew: computed beta under iteration with error precision=tolerate
% factor: estimated factor
% lambda: estimated loadings
% V: the eigenvalues matrix
% e: estimated residuals
% niter: number of interations to achieve convergence 

function [betanew, factor, lambda, V, e, niter]=Mul_betaIterNew(X, xxinv, Y, F,L, r, tolerate);
   [T,N,p]=size(X);
   betanorm=1;
   betaold=zeros(p,1);

n=0;   % number of iterations needed to stop
while (betanorm  > tolerate & n < 500) 
    n=n+1;
    [beta]=Mul_panelbetaNew(X,xxinv, Y, F,L);
    betanorm=norm(beta-betaold);
    betaold=beta;
    U=Y;
    for k=1:p;
       U=U-X(:,:,k)* beta(k);
    end

    [F,L, VNT]=panelFactorNew(U,r);
  
    
end
betanew=beta;
niter=n;
factor=F;
lambda=L;
V=VNT;
e=U-F*L'; 