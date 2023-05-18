function [cov] = matern_approx(x1,x2,Phi,lambda,Var,Lam,kernel)

%
% INPUTS
% x1 = first input values, size [N1,1]
% x2 = second input values, size [N2,1]
% Phi = eigenfunctions, size [M,1]
% lambda = eigenvalues, size [M,1]
% Var = variance parameter, size [1,1]
% Lam = 1/length-scale (shape) parameter, size [1,1]
% kernel = type of kernel
%
% 
% OUTPUTS
% cov = approximate mat covariance matrix, size [N1,N2] 

cov=zeros(length(x1),length(x2));
k = spectral_density_matern(sqrt(lambda),Var,Lam,kernel);
cov = Phi *diag(k) *Phi'; 
  

end
