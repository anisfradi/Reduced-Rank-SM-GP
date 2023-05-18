

 nfft=512; % the number of frequencies you want the spectrum evaluated at
 nperseg=50; % the number of samples to use for each spectral estimate calculation
 noverlap=40; % the number of samples to include from the calculation of spectrum N-1 in spectrum N


 gradients = {{'Black-red-yellow, custom control points', 'multigradient([0 0 0; 1 0 0; 1 1 0])'}};

 figure('DefaultAxesFontSize',18);
 tiledlayout(1,3);
 
 nexttile
 spectrogram(ynoisy,nperseg,noverlap,nfft,fs,'yaxis'); 
 gradient = gradients{1};
 colormap(gca, eval(gradient{2}));
 %colormap(parula);
 colorbar off
 title('noisy')

 
 nexttile
 spectrogram(yclean,nperseg,noverlap,nfft,fs,'yaxis'); 
 gradient = gradients{1};
 colormap(gca, eval(gradient{2}));
 %colormap(parula); 
 colorbar off
 title('clean')


 nexttile
 spectrogram(ybar,nperseg,noverlap,nfft,fs,'yaxis'); 
 gradient = gradients{1};
 colormap(gca, eval(gradient{2}));
 colorbar off
 %colormap(parula); 
 title('reconstructed')
 
 
 
 
 
