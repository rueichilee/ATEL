
% for computing the table for grand mean, time-invariant regressors, and
% common regressors. 

function table(X);


[niter, m,kk]=size(X);

TT=[5,10,15, 20, 25,        15, 20, 25];  % sample size T
NN=[100,100,100,100, 100,  15, 20,  25 ];  % number of cross-sections N 
%TT=[10,20,50, 100,100,100,100,100,100];
%NN=[100,100,100,100, 10,20,50, 100, 100];
out=zeros(kk, 12);

format=repmat( '%6.3f & ',1,12);


for i=1:kk
    a=mean(X(:,:,i));
    b=std(X(:,:,i));
    N=NN(i);
    T=TT(i);
    out(i,:)=[ a(6), b(6), a(7), b(7), a(8), b(8),a(9), b(9), a(10), b(10),a(10), b(10)] ;  % infeasible estimator
    
    fprintf(1, ['%5.0f & %5.0f &',  format, '\\\\', '\n'], N, T, out(i,:) );
end

for i=1:kk
    a=mean(X(:,:,i));
    b=std(X(:,:,i));
    N=NN(i);
    T=TT(i);
    out(i,:)=[ a(1), b(1), a(2), b(2), a(3), b(3),a(4), b(4), a(5), b(5), a(6), b(6)] ;  % interactive effect estimator
   
    fprintf(1, ['%5.0f & %5.0f &',  format, '\\\\', '\n'], N, T, out(i,:) );
end

%fid = fopen('beta_serial.table','w');
%for i=1:kk;
%    fprintf(fid, ['%5.0f & %5.0f &',  format, '\\\\', '\n'], N, T, out(i,:) );
%end
%fclose(fid)
    