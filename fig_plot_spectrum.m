

   numFreq = 200;
   ovLp = floor(numFreq/10);
   
   [pg_noisy,varpg_noisy] = welchMethod(ynoisy,numFreq,ovLp);   % computes the Periodogram (power spectrum) of 
                                                                      % a signal 'ynoisy' by 1) splitting it up
                                                                      % into sections of length 'numFreq' and with overlap 'ovLp' 2)
                                                                      % computing the fft of each section and 3) then averaging the result
                                                                      
   pg_noisy = pg_noisy/(1/2/numFreq); % convert to power spectral density 

   
   [pg_clean,varpg_clean] = welchMethod(yclean,numFreq,ovLp);  
   pg_clean = pg_clean/(1/2/numFreq);  
   
   [pg_bar,varpg_bar] = welchMethod(ybar,numFreq,ovLp);  
   pg_bar = pg_bar/(1/2/numFreq);   
   
   
  Legend={'noisy','clean','reconstructed'};

  figure('DefaultAxesFontSize',18);  
  plot(linspace(0,fs/2,length(pg_noisy)),pg_noisy,'g');
  hold on;
  plot(linspace(0,fs/2,length(pg_clean)),pg_clean,'k');
  hold on;   
  plot(linspace(0,fs/2,length(pg_bar)),pg_bar,'b');
  xlim([0 fs/2])
  set(gca,'yscale','log')
  box('off')
  xlabel('frequency (Hz)')
  ylabel('log-filter response (dB)')
  legend(Legend)
  title('power spectrum')
  