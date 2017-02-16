xdMat=zeros(4656,2);
ydMat=zeros(4656,2);
for k=2:20
    idx=force_matrix(:,k*2)~=0;
    force_temp=force_matrix(idx,k*2-1:k*2);
    disp_temp=disp_matrix(idx,k*2-1:k*2);
    no_of_points=round(length(force_temp)-(1/freq_map(k,1))*4*sampling_rate);
    xdMat(:,k)=disp_temp(no_of_points:end,2);
    ydMat(:,k)=force_temp(no_of_points:end,2);
    figure;
end