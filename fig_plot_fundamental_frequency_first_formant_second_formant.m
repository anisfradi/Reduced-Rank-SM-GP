

 % create manually the confidence interval for each formant
 I0=[50,300]; I1=[500,1000]; I2=[1500,2100];
 
 
 % find the possible values at each interval
 P0=omFit(I0(1) < omFit &  omFit < I0(2));
 P1=omFit(I1(1) < omFit &  omFit < I1(2));
 P2=omFit(I2(1) < omFit &  omFit < I2(2));
 
 
 % recover the index of possible values
 C0 = ismember(omFit, P0);
 Index0 = find(C0);
 
 C1 = ismember(omFit, P1);
 Index1 = find(C1); 
 
 C2 = ismember(omFit, P2);
 Index2 = find(C2);
 
 
 % find the fundamental frequency F0 and the two first formants F1,F2 having the maximum variances at each interval
 
 % fundamental frequency
 if length(P0)==0
     F0=NaN;
 elseif length(P0)==1
     F0=P0;
 else 
    F0=omFit(VarFit==max(VarFit(Index0)));  
 end
 
 % first formant
 if length(P1)==0
     F1=NaN;
 elseif length(P1)==1
     F1=P1;
 else 
    F1=omFit(VarFit==max(VarFit(Index1)));  
 end
 
% second formant
 if length(P2)==0
     F2=NaN;
 elseif length(P2)==1
     F2=P2;
 else 
    F2=omFit(VarFit==max(VarFit(Index2)));  
 end
 
 
 % plot the fundamental frequency F0 and the two first formants F1,F2
 figure('DefaultAxesFontSize',18);  
 stem([1 2 3],[F0,F1,F2],'b','LineWidth',3) % estimated formant values
 hold on;
 stem([1 2 3],[105,770,1810],'r--','LineWidth',2) % Praat formant values found by Praat
 grid on;
 grid minor;
 xlim([0.5 3.5])
 xticklabels({'F0','F1','F2'})
 ylabel('frequency (Hz)')
 legend('estimated','praat')
 title('fundamental frequency and two first formant values')

 
 
 


