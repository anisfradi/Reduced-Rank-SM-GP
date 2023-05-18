function [varx,lamx,omx,vary,Info] = fit_SM_GP_hyperparameter_temporal_domain(x,y,Lambda,Phi,D,kernel,fs,varargin)


%
% INPUTS
% x = time instances, size [N,1]
% y = the noisy signal, size [N,1]
% Lambda = eigenvalues, size [M,1]
% Phi = eigenfunctions, size [M,1]
% D = number of quasi-periodic components
% kernel = type of kernel
% fs = sampling frequency
% OPTIONAL INPUTS:
% Opts = set of options
% OUTPUTS
% varx = variance parameter, size [D,1] 
% lamx = 1/length-scale (shape) parameter, size [D,1] 
% om = center filter bank parameter, size [D,1]
% vary = noise variance, size [1,1]
% Info = structure containing information about the estimation
%        process 



% Rescale and zero-mean y (we rescale the variances at the end to compensate)
x = x(:); % Ensure column vectors and center data 
y = y(:)- mean(y(:));
varSig = var(y);
y = y/sqrt(varSig);
 
% parameter initialization
Var=repmat(var(y)/D,D,1);
Lam=repmat(100,D,1);
om=linspace(10,fs/2-fs/10,D)';
vary=log(0.1);

theta=[log(Var)  log(Lam)  log(om)];
theta=[theta(:)  ; vary];
  
numIts = varargin{1}.numIts;
% Learning hyperparameters
[theta,Obj,inCur] = minimize(theta,'get_Obj_matern_temporal',numIts, ...
				  x,y,Lambda,Phi,D,kernel);


varx=exp(theta(1:D));
lamx=exp(theta(D+1:2*D));
omx=exp(theta(2*D+1:3*D));  
vary=exp(theta(end));


Info.Objs = Obj;
Info.ins = inCur;


end

