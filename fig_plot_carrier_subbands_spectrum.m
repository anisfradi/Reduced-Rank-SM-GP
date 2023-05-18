

numFreq = N;
ovLp = 0;
      
[pg_clean,varpg_clean] = welchMethod(yclean,numFreq,ovLp);  
pg_clean = pg_clean/(1/2/numFreq); 
   
   
for d=1:D
    spectrum = welchMethod(fbar(:,d),numFreq,ovLp);  
    fbar_spectrum(:,d) = spectrum /(1/2/numFreq);
end

[time,frequency]=meshgrid(linspace(0,N/fs,N),omFit);


figure('DefaultAxesFontSize',18);
 
tiledlayout(2,1);

nexttile
p1=waterfall(time,frequency,fbar_spectrum');
%set(p1, 'EdgeColor', [0 0.4470 0.7410]);
set(gca,'ZScale','log');
ylabel('frequency (Hz)')
grid on
grid minor
xticks([])
%yticks([])
zticks([])
xlim([0,N/fs])
ylim([0,fs/2])
title('log-spectrum GP carrier subbands')
view(5,80)



nexttile
plot(linspace(0,fs/2,length(pg_clean)),pg_clean,'k')
xlabel('frequency')
xticks([])
yticks([])
xlim([0,fs/2])
set(gca,'yscale','log')
title('log-spectrum clean signal')

