

clear;
close all;
clc;


%% Load signal, downsample if fs > 8 Khz and normalize it to unit variance
File = 'had3-nat'; % Name of file to load
fs_ = 8000; % sampling rate of file
y_max_length = 1; % seconds
D = 20; % number of quasi-periodic components
[x,yclean,fs,N] = from_audio_to_data(File,y_max_length,fs_); % gives: x=inupts, yclean=clean signal,
                                                             % fs=frequency sampling and N=sample size

%% add a noise to the signal 'yclean'
vary=0.1;  % fix the noise variance
ynoisy = get_noisy_signal(yclean,N,vary); % noisy signal


%% kernel type
kernel='exp'; 


%% hyperparameter learning in the temporal domain by Fradi (optional)                                                                                                           
M=100; % truncation order of RR-GP
[eigenf,eigenv] = eigen_val_fct(x,M); % evaluate eigenfunctions and eigenvalues of the Laplace operator
opts.numIts = 50; % number of iterations needed to minimize the objective function
[VarFit,LamFit,omFit,varyFit,Info] = fit_SM_GP_hyperparameter_temporal_domain(x,ynoisy,eigenv,eigenf,D,kernel,fs,opts); % find the optimal hyperparameter minimizing
                                                                                                                % the minus marginal likelihood in the spectral domain                                                                                                                       % the minus marginal likelihood in the temporal domain
                                                                                                        
%% plot the fitted variance parameter VarFit vs the center parameter omFit
fig_plot_variance_center_parameter


%% estimate the fundamental frequency F0, the first formant F1 and second formant F2 and compare 
%% them  to values given by praat 
fig_plot_fundamental_frequency_first_formant_second_formant


%% compute the posterior mean of each component
M=100; % truncation order of RR-GP
[eigenf,eigenv] = eigen_val_fct(x,M); % evaluate eigenfunctions and eigenvalues of the Laplace operator
ynoisy_centered=ynoisy-mean(ynoisy);  % remove the mean to get a zero-mean signal
[fbar,covbar]=posterior_mean_covariance_matern_approx(x,ynoisy_centered,eigenf,eigenv,VarFit,LamFit,omFit,kernel,varyFit,D); % using the Cholesky decomposition
                                                                                                                             % for an efficient implementation and
                                                                                                                             % numerical stability
%% plot the GP carrier subbands                                                                                                                
fig_plot_carrier_subbands                                                                                                                


%% plot the log-spectrum GP carrier subbands                                                                                                                
fig_plot_carrier_subbands_spectrum    


%% compute the posterior mean of the signal
ybar=sum(fbar,2)+mean(ynoisy); % add the mean removed before to get a non zero-mean reconstructed signal


%% plot the power log-spectrum 
fig_plot_spectrum


%% plot the spectrogram 
fig_plot_spectrogram

