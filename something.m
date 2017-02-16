%% Adjusting with the log file frequency data
TT=VarName11;
[Y, M, D, H, MN, S] = datevec(TT);
TT=H*3600+MN*60+S;
log=[VarName3,TT];
freq_comp=freq_detected;
freq_adj=zeros(length(freq_comp),2);
freq_adj(:,2)=freq_comp(:,2);
for j=3:length(log)
    temp_idx=freq_comp(:,2)<log(j,2)&freq_comp(:,2)>log(j-1,2);
    freq_adj(temp_idx,1)=log(j-1,1);
end
%% 
check=round(smooth(freq_detected(:,1),20),2);
check2=round(smooth(check,20),2);
plot(freq_detected(:,2),check2,'+')
hold on
plot(freq_detected(:,2),check1,'o')