function [opts] = option_list()

  opts.verbose = 1; % view plots of the fitting process
  opts.minT = 100;
  opts.maxT = 1000;
  opts.numIts = 10; % number of iterations
  opts.numLevels = 30; % number of levels
  opts.bet = 750;  % increase if filters bunch together
  opts.reassign = 0;  % set to 1 if filters really bunch together

end

