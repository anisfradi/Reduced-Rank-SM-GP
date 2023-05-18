function [S1,S2] = spectral_density_matern_gradient(w,Var,Lam,kernel)


%
% INPUTS
% w = frequency value
% Var = variance parameter, size [1,1]
% Lam = 1/length-scale (shape) parameter, size [1,1]
% kernel = kernel type 
%
% 
% OUTPUTS
% S1 = partial derivative w.r.t variance parameter of spectral density of matérn, size [1,1] 
% S2 = partial derivative w.r.t 1/length-scale (shape) parameter of spectral density of matérn, size [1,1] 



        if strcmp(kernel,'exp')
           Lam = Lam;
           nu=1/2;
           S2=2*Var*sqrt(pi)*(gamma(nu + 1/2)/gamma(nu))*((2*nu*Lam^(2*nu-1))./(Lam * Lam + w.^2).^(nu + 1/2) -...
       (2*(nu+1/2)*Lam^(2*nu+1))./(Lam * Lam + w.^2).^(nu +3/2));      
        elseif strcmp(kernel,'matern32')
           Lam = sqrt(3) .* Lam;
           nu=3/2;
           S2=2*Var*sqrt(pi)*(gamma(nu + 1/2)/gamma(nu))*((2*nu*Lam^(2*nu-1))./(Lam * Lam + w.^2).^(nu + 1/2) -...
       (2*(nu+1/2)*Lam^(2*nu+1))./(Lam * Lam + w.^2).^(nu +3/2));     
        elseif strcmp(kernel,'matern52')
           Lam = sqrt(5) .* Lam;
           nu=5/2;
           S2=2*Var*sqrt(pi)*(gamma(nu + 1/2)/gamma(nu))*((2*nu*Lam^(2*nu-1))./(Lam * Lam + w.^2).^(nu + 1/2) -...
       (2*(nu+1/2)*Lam^(2*nu+1))./(Lam * Lam + w.^2).^(nu +3/2));  
        elseif strcmp(kernel,'rbf')
           S2=2*Var*sqrt(pi)*(-1./(Lam*Lam) + (2*w.^2)./(Lam*Lam*Lam*Lam)).*exp(-w.^2./(Lam*Lam));
        end
      
        
       S1=spectral_density_matern(w,Var,Lam,kernel)/Var;      

end

