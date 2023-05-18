function [fbar,cov] = posterior_mean_covariance_matern_approx(x,y,Phi,lambda,Var,Lam,om,kernel,var_y,D)

%
% INPUTS
% x = time instances, size [N,1]
% y = observed signal, size [N,1]
% Phi = eigenfunctions, size [M,1]
% lambda = eigenvalues, size [M,1]
% Var = variance parameter, size [D,1]
% Lam = 1/length-scale (shape) parameter, size [D,1]
% om = center filter bank parameter, size [D,1]
% kernel = type of kernel
% var_y = noise variance 
% D = number of quasi-periodic components
%
% 
% OUTPUTS
% fbar = approximate posterior mean of D components with mixture of approx mat, size [N,D] 
% cov = approximate posterior covariance of of D components with mixture of approx mat, size [N,N,D] 

N=length(x);

C_train=Mixture_matern_approx(x,x,Phi,lambda,Var,Lam,om,kernel,D);

% Calculate the Cholesky factor for fast efficient implementation and numerical stability 
L=chol(C_train+var_y*eye(N),'lower');  % lower triangular matrix
Sigma_y_inv=pinv(L')*pinv(L); % the inversed signal covariance matrix 
  
  
fbar=[];   

cov=[];

for d=1:D
    A=cos(om(d)*(x'-x)).*matern_approx(x,x,Phi,lambda,Var(d),Lam(d),kernel); % prior covariance of d-th component
    fbar(:,d)=A*Sigma_y_inv*y; % d-th posterior mean 
    cov(:,:,d)= var_y*A*Sigma_y_inv; % d-th covariance matrix
end


end