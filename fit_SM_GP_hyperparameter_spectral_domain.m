function [varx,lamx,om,vary,Info] = fit_SM_GP_hyperparameter_spectral_domain(y,D,kernel,fs,vary,varargin)


%
% INPUTS
% y = the noisy signal, size [N,1]
% D = number of quasi-periodic components
% kernel = type of kernel
% fs = sampling frequency
% vary = noise variance initialization
% OPTIONAL INPUTS:
% Opts = set of options
% 
% OUTPUTS
% varx = variance parameter, size [D,1] 
% lamx = 1/length-scale (shape) parameter, size [D,1] 
% om = center filter bank parameter, size [D,1]
% vary = noise variance, size [1,1]
% Info = structure containing information about the estimation
%        process including
%        Objs = objectives at each scale concatenated into a vector
%        ins = number of completed iterations at each scale


% Rescale and zero-mean y (we rescale the variances at the end to
% compensate)
y = y(:) - mean(y(:));
varSig = var(y);
y = y/sqrt(varSig);




if ~isfield(varargin{1},'theta_init')
  % Initialise the STFT parameters
    % uniformly over log-frequency and make them fairly broad
    mVar = ones(D,1)/D;
    FLim = [1/40,0.35]; % limits for the centre frequency initialisation
    dfFrac = 1/20; % width of the processes wrt centre-frequency
%     fmax = logspace(log10(FLim(1)),log10(FLim(2)),D)';
    fmax = linspace(FLim(1),FLim(2),D)';
    [om,lamx] = freq2probSpec(fmax,fmax*dfFrac,1);

  cvar_d = mVar .* (1 - lamx.^2); % conditional variance
  % convert hypers to continuous form so we can calculate spectral density
  lam_c = -log(lamx);
  lam_max = 0.4;
  lamLim = ones(D,1)*[0,lam_max];
%   lam_c = min(-log(lamx), lam_max-1e-5);
  cvar_c = cvar_d ./ (1 - exp(-2 .* lam_c));
  mVar_c = cvar_c ./ (1 - lam_c.^2);
%   omLim = ones(D,1)*[0.0,pi]; % limits on the centre frequency
else
  theta_init = varargin{1}.theta_init;
  cvar_c = theta_init(1:D);
  if strcmp(kernel,'exp')
    lam_c = theta_init(D+1:2*D);
  elseif strcmp(kernel,'matern32')
    lam_c = theta_init(D+1:2*D) .* sqrt(3);
  elseif strcmp(kernel,'matern52')
    lam_c = theta_init(D+1:2*D) .* sqrt(5);
  end
  om = theta_init(2*D+1:end);
%   cvar_d = mVar .* (1 - lamx.^2); % conditional variance
  % convert hypers to continuous form so we can calculate spectral density
  if ~isfield(varargin{1},'bandwidth_lim')
      bandwidth_lim = 2;
  else
      bandwidth_lim = varargin{1}.bandwidth_lim;
  end
  lam_max = min(lam_c .* bandwidth_lim,1-1e-5);
%   lam_c = min(-log(lamx), lam_max-1e-5);
%   cvar_c = cvar_d ./ (1 - exp(-2 .* lam_c));
  mVar_c = max(cvar_c ./ (1 - lam_c.^2),1e-3);
  lamLim = [zeros(D,1),lam_max];
%   omLim = [om-0.1,om+0.1]; % limits on the centre frequency
  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constraints on the variables
minVar = max(mVar_c/400,1e-5); % minimum marginal variance - might want to taper
omLim = ones(D,1)*[0.0,pi]; % limits on the centre frequency

% if strcmp(kernel,'matern52')
%     keyboard
%   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COARSE TO FINE PROCESS

if ~isfield(varargin{1},'numLevels')
  numLevels = 40; % number of levels
else
  numLevels = varargin{1}.numLevels;
end

if ~isfield(varargin{1},'numIts')
  numIts = 10; % number of iterations per level
else
  numIts = varargin{1}.numIts;
end

if ~isfield(varargin{1},'minT')
  minT = 200;%200; % minimum segment size to compute spectrum of
else
  minT = varargin{1}.minT;
end

if ~isfield(varargin{1},'maxT')
  maxT = 1000; % maximum segment size to compute spectrum of
else
  maxT = varargin{1}.maxT;
end




% Observation noise tapered
 if ~isfield(varargin{1},'vary_an')
   vary_an = logspace(log10(1e-6),log10(1e-10),numLevels); % annealing settings
 else
   vary_an = varargin{1}.vary_an;
 end


if ~isfield(varargin{1},'bet')
  bet = logspace(log10(100),0,numLevels); % how strongly to encourage variances to be pruned
else
  bet = logspace(log10(varargin{1}.bet),0,numLevels);
%  bet = varargin{1}.bet;
end


		 
T = length(y);

numFreq = floor(logspace(log10(min([minT,T])),log10(min([maxT,T])),numLevels));
ovLp = floor(numFreq/10);

Objs = [];

% keep the historical values of the variables as they evolve
omHist = repmat(NaN,[D,numLevels+1]);
lamHist = repmat(NaN,[D,numLevels+1]);
mVarHist = repmat(NaN,[D,numLevels+1]);

omHist(:,1) = om;
lamHist(:,1) = lam_c;
mVarHist(:,1) = mVar_c;

% Compute the likelihood if the user has requested it
if nargin>3
  if isfield(varargin{1},'yHO')
    compHOLike = 1;
    
    % compute the full periodogram for the held out data
    yHO = varargin{1}.yHO;
    yHO = yHO(:)/sqrt(varSig);
    
    THO = length(yHO);
    likeHO = repmat(NaN,[numLevels,1]);
	
    [pgHO,varpg] = welchMethod(yHO,THO,0);
    pgHO = pgHO/(1/2/THO);
    
    [pgUR,varpg] = welchMethod(y,THO,0);
  
    if mod(THO,2)==0
      % if even
      specHO = [pgHO;pgHO(end-1:-1:2)];
    else
      % if odd
      specHO = [pgHO;pgHO(end:-1:2)];
    end
  else
    compHOLike = 0;
  end
else
  compHOLike = 0;
end

 
% compute the full periodogram for the unregularised objective
likeUnReg = repmat(NaN,[numLevels,1]);

[pgUR,varpg] = welchMethod(y,T,0);
pgUR = pgUR/(1/2/T);
    
if mod(T,2)==0
  % if even
  specUR = [pgUR;pgUR(end-1:-1:2)];
else
  % if odd
  specUR = [pgUR;pgUR(end:-1:2)];
end
    




for c2f=1:numLevels

  % Estmate spectrum of y using Welch's periodogram
  [pg,varpg] = welchMethod(y,numFreq(c2f),ovLp(c2f));
  
  
  pg = pg/(1/2/numFreq(c2f)); % convert to power spectral density estimate
  

  %pg = pg*numFreq;

  if mod(T,2)==0
    % if even
    specTar = [pg;pg(end-1:-1:2)];
  else
    % if odd
    specTar = [pg;pg(end:-1:2)];
  end  
  
%   minVar = var_c ./ 2.;
%   minVar(:) = 0.001;
  
  % Constraints on variables
%   vary = vary_an(c2f);

  % current settings of the variables
  theta = [log(mVar_c-minVar);...
	   log(om-omLim(:,1))-log(omLim(:,2)-om);...
	   log(lam_c-lamLim(:,1))-log(lamLim(:,2)-lam_c);...
       log(vary)]; 
  % Base line noise to avoid divide by zero problems
  
  Obj_func = strcat('get_Obj_',kernel,'_spectral')
  Obj_func_handle = str2func(Obj_func);
  
  
name = 'Level %d\n'; fprintf(name,c2f); % print the level
  
  
  % Fit the spectrum using probabilistc spectrogram
  [theta,ObjCur,inCur] = minimize(theta,Obj_func,numIts, ...
				  specTar,minVar,omLim, ...
				  lamLim,bet(c2f)*numFreq(c2f)/numFreq(c2f(1)));


  if compHOLike==1
    % Compute HO likelihood if asked for
    [likeHO(c2f,1),dObjTemp] = Obj_func_handle(theta,specHO,minVar, ...
						   omLim,lamLim,0);    
  end
  
  % Compute the real objective on the full data too without regularisation
  [likeUnReg(c2f,1),dObjTemp] = Obj_func_handle(theta,specUR,minVar, ...
							 omLim,lamLim,0);
  
  
  % Collect information about the current iteration
  Objs = [Objs;ObjCur];
  ins(c2f) = inCur;
  
  % Compute the current settings of the parameters
  mVar_c= minVar+exp(theta(1:D));
  om = omLim(:,1)+(omLim(:,2)-omLim(:,1))./(1+exp(-theta(D+1:2*D)));
  lam_c = lamLim(:,1)+(lamLim(:,2)-lamLim(:,1))./(1+exp(-theta(2*D+1:3*D)));
  vary = exp(theta(end));
  
  
  vary = max(vary_an(c2f),vary); % At each level the noise variance takes the maximum between the observation noise tapered 
                                 % and the last estimated noise variance to avoid divisions by zero
    
  
  % Store historical values for the processes
  omHist(:,c2f+1) = om;
  lamHist(:,c2f+1) = lam_c;
  mVarHist(:,c2f+1) = mVar_c;
  mVarHist(:,c2f+1) = vary;

  % plot if requested by user
  
  if nargin>3
    if isfield(varargin{1},'verbose')
      if varargin{1}.verbose==1
	% plot the fit
	if compHOLike==1
	  figH1 =  plot_fit_probSTFT_kern_cts(pg,mVar_c,om,lam_c,mVarHist,omHist,lamHist, ...
				    Objs,ObjCur,minVar,kernel,likeHO,likeUnReg);  
	else
	  figH1 =  plot_fit_probSTFT_kern_cts(pg,mVar_c,om,lam_c,mVarHist,omHist,lamHist, ...
				    Objs,ObjCur,minVar,kernel);  
	end
      end
    end
  end

% % re-initialised the pruned components of the model - don't reassign
% on the last iteration

if isfield(varargin{1},'reassign')&c2f<numLevels-1
  
  if ~isfield(varargin{1},'reassignThresh')
    reassignThresh = 8;
  else
    reassignThresh = varargin{1}.reassignThresh;
  end
  
  if varargin{1}.reassign==1
    for d=1:D    
      if mVar_c(d)/minVar(d)<reassignThresh

	disp('moving pruned process')
	  varx = mVar_c.*(1-lam_c.^2);
	  
	  % remove current process
	  lamxCur = lam_c; lamxCur(d) = [];
	  varxCur = varx; varxCur(d) = [];
	  omCur = om; omCur(d) = [];

	  N = length(specTar);
	  freqs = linspace(0,0.5,ceil(N/2));
	
      specModHandle = str2func(strcat('get_pSTFT_spec_cts_',kernel));
	  specMod = sum(specModHandle(freqs,lamxCur,varxCur,omCur));
	  
	  
	  % reassigning based on differences in log-spectra
	  dspec = log(specTar(1:ceil(N/2)))-log(specMod');
	    
	  [val,pos] = max(dspec);
	  
	  mVar_c(d) = 1/20;
	  om(d) = 2*pi*freqs(pos);
	  lam_c(d) = 0.05;
	  varx = mVar_c.*(1-lam_c.^2);
	  
      
      end
    end
  end
end



end

Info.Objs = Objs;
Info.ins = ins;
Info.dObj=dObjTemp;


% we recover the conditional variances based on the marginal variances so that
% their sum is equal to that of the signal 

lamx = lam_c;
mVar = mVar_c;
rescale = varSig/sum(mVar);
varx = mVar.*(1-lamx.^2);
varx = rescale*varx;

om=om*fs/(2*pi); % from [0,pi] to [0,fs/2]
lamx=lamx*fs/2; % from [0,1] to [0,fs/2]

if compHOLike==1
  Info.likeHO=likeHO;
end

Info.likeUnReg=likeUnReg;

