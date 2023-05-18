function [y] = get_noisy_signal(y,N,var)

y=y+mvnrnd(repmat(0,N,1),diag(repmat(var,N,1)))'; 

end

