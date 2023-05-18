function [eigenf,eigenv] = eigen_val_fct(x,M)

Lt = 1;   % Boundary, x \in [-Lt,Lt]
NN = (1:M)'; % Indices for eigenvalues
 
%% Define eigenvalue and eigenfunction and of the Laplacian in 1D
eigenval = @(n) (n(:)'*pi/2/Lt).^2;
eigenfun = @(n,x) Lt^(-1/2)*sin(kron(n(:)'*pi,(x(:)+Lt)/2/Lt));        
%%  EValuate eigenfunctions and eigenvalues       
eigenf = eigenfun(NN,x); 
eigenv = eigenval(NN)';  


end

