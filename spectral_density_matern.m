function [S] = spectral_density_matern(w,Var,Lam,kernel)


%
% INPUTS
% w = frequency value
% Var = variance parameter, size [1,1]
% Lam = 1/length-scale (shape) parameter, size [1,1]
% kernel = kernel type 
%
% 
% OUTPUTS
% S = spectral density of matérn, size [1,1] 


        if strcmp(kernel,'exp')
           Lam = Lam;
           nu=1/2;
           S=(2*Var*sqrt(pi)*Lam^(2*nu)*gamma(nu + 1/2)/gamma(nu))./(Lam * Lam + w.^2).^(nu + 1/2);
        elseif strcmp(kernel,'matern32')
           Lam = sqrt(3) .* Lam;
           nu=3/2;
           S=(2*Var*sqrt(pi)*Lam^(2*nu)*gamma(nu + 1/2)/gamma(nu))./(Lam * Lam + w.^2).^(nu + 1/2);    
        elseif strcmp(kernel,'matern52')
           Lam = sqrt(5) .* Lam;
           nu=5/2;
           S=(2*Var*sqrt(pi)*Lam^(2*nu)*gamma(nu + 1/2)/gamma(nu))./(Lam * Lam + w.^2).^(nu + 1/2);
        elseif strcmp(kernel,'rbf')
           S=2*Var*sqrt(pi/(Lam * Lam))*exp(-w.^2/(Lam * Lam));
        end
                     

end

