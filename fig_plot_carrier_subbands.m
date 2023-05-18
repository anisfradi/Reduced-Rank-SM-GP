
[time,frequency]=meshgrid(linspace(0,N/fs,N),omFit);


figure('DefaultAxesFontSize',18);
 
tiledlayout(2,1);

nexttile
p1=waterfall(time,frequency,fbar');
%set(p1, 'EdgeColor', [0 0.4470 0.7410]);
ylabel('frequency (Hz)')
grid on
grid minor
xticks([])
%yticks([])
zticks([])
xlim([0,N/fs])
ylim([0,fs/2])
title('GP carrier subbands')
view(5,80)

nexttile
plot(linspace(0,N/fs,N),yclean,'k')
xlabel('time')
xticks([])
yticks([])
xlim([0,N/fs])
title('clean signal')