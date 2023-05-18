function [cov] = Mixture_matern_approx(x1,x2,Phi,lambda,Var,Lam,om,kernel,D)


%
% INPUTS
% x1 = first input values, size [N1,1]
% x2 = second input values, size [N2,1]
% Phi = eigenfunctions, size [M,1]
% lambda = eigenvalues, size [M,1]
% Var = variance parameter, size [D,1]
% Lam = 1/length-scale (shape) parameter, size [D,1]
% om = center filter bank parameter, size [D,1]
% kernel = type of kernel
% D = number of quasi-periodic components
%
% 
% OUTPUTS
% cov = mixture of approximate mat covariance matrix, size [N1,N2] 


cov=zeros(length(x1),length(x2));

for n=1:D
    
    cov=cov+cos(om(n)*(x1'-x2)).*matern_approx(x1,x2,Phi,lambda,Var(n),Lam(n),kernel);

end

end
