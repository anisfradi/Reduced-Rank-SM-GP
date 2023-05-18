function [Obj,varargout] = get_Obj_exp_spectral(theta,specTar,minVar,limOm,limLam,bet)


%
% INPUTS
% theta = vector of parameters of be estimated, size [3*D+1,1]
% containing: variance parameter, size [D,1]; 1/length-scale (shape) parameter, size [D,1]; center filter bank parameter, size [D,1]; noise variance, size [1,1]


% specTar = target spectrum of the data to be fit, size [N,1]
% minVar = minimum of the marginal variances, size [D,1]
% limOm = min/max centre-frequencies e.g. [0,1/2], size [D,2]
% limLam = min/max bandwidths e.g. [0,1], size [D,2]
% bet = strength of the gamma shrinkage prior on the marginal
%       variance parameter (set to zero for no pruning)
% 
% OUTPUTS
% Obj = objective
% OPTIONAL OUTPUT:
% dObj = derivative of the objective wrt the parameters, size [3*D+1,1]
% 


D = (length(theta)-1)/3;

dVar= exp(theta(1:D));
mVar = minVar+dVar;

om = limOm(:,1)+(limOm(:,2)-limOm(:,1))./(1+exp(-theta(D+1:2*D)));
lam = limLam(:,1)+(limLam(:,2)-limLam(:,1))./(1+exp(-theta(2*D+1:3*D)));
vary = exp(theta(end));
% cts model parameters:


N = length(specTar);
omegas = linspace(0,pi,ceil(N/2));
omegas = [omegas,-omegas([floor(N/2):-1:1])];

% conditional variance
cVar = mVar .* (1 - lam.^2);

% Get the component spectra
spec = ones(1,N)*vary;
for d=1:D

  % for the objective
  alp1_c = lam(d).^2 + (omegas-om(d)).^2;
  alp2_c = lam(d).^2 + (omegas+om(d)).^2;
  spec = spec + cVar(d) .* lam(d) .* (alp1_c.^-1 + alp2_c.^-1);
    
end


% Minus of marginal log-likelihood in the spectral domain (cost function)
Obj1 = sum(log(spec))+sum(specTar'./spec);

% Prior cost function (term of regularization)
Obj2 = bet*sum(mVar);

% Rescale the regularized cost function to get nice units (objectif function to be minimized)
Obj = (Obj1+Obj2)/N;

%plot(log(specTar),'-k'); hold on; plot(log(spec),'-r')
%keyboard

if nargout>1
  % Derivative of the components wrt parameters initialized to zeros

  dObjdtransVar = zeros(D,1);
  dObjdtransOm = zeros(D,1);
  dObjdtransLam = zeros(D,1);
  
  dObjdspec = 1./spec-specTar'./(spec.^2);  % this is a common term yielding for all partial derivatives to be multiplied by the derivative of spec wrt each parameter
  
 
  % A practical consideration is that hyperparameter optimization benefits from parametrising the model with the marginal variances rather than the conditional variances
  for d=1:D
 
    % Derivative wrt transformed marginal variance
    alp1_c = lam(d).^2 + (omegas-om(d)).^2;
    alp2_c = lam(d).^2 + (omegas+om(d)).^2; 
    dspecdtransVar = (mVar(d)-minVar(d)) .* (1 - lam(d).^2) .* lam(d) .* (alp1_c.^-1 + alp2_c.^-1);
    dObjdtransVar(d) = sum(dObjdspec.*dspecdtransVar); 
    
    
    % Derivative wrt transformed centre frequency
    dspecdom = 2 .* mVar(d) .* (1 - lam(d).^2) .* lam(d) * (alp1_c.^-2 .* (omegas-om(d)) - alp2_c.^-2 .* (omegas+om(d)));
    domdtransOm = (limOm(d,2)-limOm(d,1))*(1/4)./cosh(theta(D+d)/2).^2;
    dObjdtransOm(d) = sum(dObjdspec.*dspecdom)*domdtransOm; 

    
    % Derivative wrt transformed lambda
%     dspecdlam = mVar(d) .* (alp1_c.^-1 + alp2_c.^-1 - 2 .* lam(d).^2 * ...
%                           (alp1_c.^-2 + alp2_c.^-2));
    dspecdlam = mVar(d) .* ((1 - 3.*lam(d).^2) .* (alp1_c.^-1 + alp2_c.^-1) ...
                            - 2 .* lam(d).^2 .* (1 - lam(d).^2) * (alp1_c.^-2 + alp2_c.^-2) );
    dlamdtransLam = (limLam(d,2)-limLam(d,1))*(1/4)./cosh(theta(2*D+d)/2).^2;
    dObjdtransLam(d) = sum(dObjdspec.*dspecdlam)*dlamdtransLam;  
  end
 
  % Derivative wrt noise variance
  
  dObjdtransvary = sum(dObjdspec)*1e-0; 
  
  
  dObj1 = [dObjdtransVar;dObjdtransOm;dObjdtransLam;dObjdtransvary];

  dObj2 = [bet*dVar;zeros(2*D+1,1)];
  
  dObj = [dObj1+dObj2]/N;
  varargout{1}= dObj;
  
end

