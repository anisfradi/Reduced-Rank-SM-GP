function [Obj, dObj] = get_Obj_matern_temporal(theta,x,y,lambda,Phi,D,kernel)


%
% INPUTS
% theta = vector of parameters of be estimated, size [3*D+1,1]
% containing: variance parameter, size [D,1]; 1/length-scale (shape) parameter, size [D,1]; center filter bank parameter, size [D,1]; noise variance, size [1,1]

% x = input values, size [N,1]
% y = the observed signal, size [N,1]
% Phi = eigenfunctions, size [M,1]
% lambda = eigenvalues, size [M,1]
% D = number of quasi-periodic components
% kernel = type of kernel

% 
% OUTPUTS
% Obj = objective
% OPTIONAL OUTPUT:
% dObj = derivative of the objective wrt the parameters, size [3*D,1]
% 


  N = numel(y);  % Number of N=observations 
  

  % Extract parameters
  Var = exp(theta(1:D)); %  variance parameter
  Lam  = exp(theta(D+1:2*D)); % 1/length-scale (shape) parameter
  omega =  exp(theta(2*D+1:3*D)); % center filter bank parameter
  sigma2 = exp(theta(end));

  
 % Initialize the covariance matrix 
  cov=zeros(N,N); 

 for d=1:D
  % Evaluate the covariance matrix
  k = spectral_density_matern(sqrt(lambda),Var(d),Lam(d),kernel);
  PhiPhi = Phi *diag(k) *Phi'; % approximate covariance matrix     
  cov=cov+cos(omega(d)*(x-x')).*PhiPhi;
 end 
 
  % Calculate the Cholesky factor
  L=chol(cov + sigma2*eye(N),'lower');  
  Kinv=pinv(L')*pinv(L);  
  A=2*sum(log(diag(L)));
  alpha = L'\(L\y); 
  B=y'*alpha;
  % Return negative log marginal likelihood in the temporal domain
  Obj=.5*(A+B) + .5*N*log(2*pi);
  
  
  % Initialize the partial derivatives of negative log marginal likelihood
  dObj_Var = zeros(D,1); 
  dObj_Lam = zeros(D,1);
  dObj_om = zeros(D,1); 
  dObj_vary=zeros(1,1);
  
  for d=1:D
      % Evaluate the partial derivatives of spectral density
      [kk1 , kk2]=spectral_density_matern_gradient(sqrt(lambda),Var(d),Lam(d),kernel);
      
      dPhiPhi_Var=Phi *diag(kk1) *Phi'; % derivative of covariance matrix w.r.t variance parameter
      dPhiPhi_Lam=Phi *diag(kk2) *Phi'; % derivative of covariance matrix w.r.t 1/length-scale (shape) parameter
    
      dcov_Var=cos(omega(d)*(x-x')).*dPhiPhi_Var;
      dcov_Lam=cos(omega(d)*(x-x')).*dPhiPhi_Lam;
      
     % Evaluate the partial derivatives of negative log marginal likelihood
      dObj_Var(d)=-.5*trace(( alpha*alpha' - Kinv) *dcov_Var);   
              
      dObj_Lam(d)=-.5*trace(( alpha*alpha' - Kinv) *dcov_Lam);   
     
       dObj_om(d)=-.5*trace(( alpha*alpha' - Kinv) *(-(x-x').*sin(omega(d)*(x-x')) .* PhiPhi));   
  
%     % Account for the log-transformed values
       dObj_Var(d) = exp(theta(d))*dObj_Var(d);      
       dObj_Lam(d) = exp(theta(D+d))*dObj_Lam(d);
       dObj_om(d) = exp(theta(2*D+d))*dObj_om(d);
 
  end
  
  d_Obj_vary=.5*trace(Kinv) -.5*y'*Kinv*Kinv*y;
  
  d_Obj_vary= exp(theta(end))*d_Obj_vary;
  
    
  dObj = [dObj_Var ; dObj_Lam ; dObj_om ; d_Obj_vary];
   
end
