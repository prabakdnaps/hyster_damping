clear all;
close all;
tic;
% %Constants
% freq_resolution=0.05;
% start_freq=42;
% end_freq=55;
% sampling_rate=3200;
% start_time=3.3;
% end_time=480;
% no_of_periods=4;
% voltage =  importdata('E:\Thesis\aorta\15cm_100mmHg_nonlinear\0.45N\0.45N_up_volt.txt','\t',31);
% force =  importdata('E:\Thesis\aorta\15cm_100mmHg_nonlinear\0.45N\0.45N_up_force.txt','\t',31);
% disp =  importdata('E:\Thesis\aorta\15cm_100mmHg_nonlinear\0.45N\0.45N_up_disp.txt','\t',31);

%Constants
% freq_resolution=0.05;
% start_freq=50;
% end_freq=60;
% sampling_rate=3200;
% start_time=4;
% end_time=555;
% no_of_periods=4;
% voltage =  importdata('E:\Thesis\aorta\15cm_100mmHg_nonlinear\0.05N\0.05N_up_volt_low.txt','\t',31);
% force =  importdata('E:\Thesis\aorta\15cm_100mmHg_nonlinear\0.05N\0.05N_up_force.txt','\t',31);
% disp =  importdata('E:\Thesis\aorta\15cm_100mmHg_nonlinear\0.05N\0.05N_up_disp.txt','\t',31);

freq_resolution=0.05;
start_freq=82;
end_freq=110;
sampling_rate=6400;
start_time=4;
end_time=970;
no_of_periods=4;
force =  importdata('E:\Thesis\aorta\test_cases\Rubber_Panel\2N_UP\2N_UP_force.txt','\t',34);
disp =  importdata('E:\Thesis\aorta\test_cases\Rubber_Panel\2N_UP\2N_UP_disp.txt','\t',34);

% freq_resolution=0.05;
% start_freq=98;
% end_freq=108;
% sampling_rate=6400;
% start_time=4;
% end_time=135;
% no_of_periods=4;
% force =  importdata('E:\Thesis\aorta\test_cases\Rubber_Panel\0.5N_UP\0.5N_RP_up_force.txt','\t',34);
% disp =  importdata('E:\Thesis\aorta\test_cases\Rubber_Panel\0.5N_UP\0.5N_RP_up_disp.txt','\t',34);


% freq_resolution=0.1;
% start_freq=15;
% end_freq=25;
% sampling_rate=1600;
% start_time=4;
% end_time=2700;
% no_of_periods=4;
% force =  importdata('E:\Thesis\aorta\test_cases\Rubber_Plate\0.01N\0.01N_up_force.txt','\t',34);
% disp =  importdata('E:\Thesis\aorta\test_cases\Rubber_Plate\0.01N\0.01N_up_disp.txt','\t',34);

voltage = force.data;
start_time_idx=voltage(:,1)>start_time & voltage(:,1)<end_time;
voltage=voltage(start_time_idx,:);
force=force.data;
temp_idx=force(:,1)<voltage(end,1)&force(:,1)>voltage(1,1);
force=force(temp_idx,:);
disp=disp.data;
temp_idx=disp(:,1)<voltage(end,1)&disp(:,1)>voltage(1,1);
disp=disp(temp_idx,:);
plot(voltage(:,1),voltage(:,2));
save('voltage.mat','voltage');
save('force.mat','force');
save('disp.mat','disp');
toc;
% plot(voltage(:,1),voltage(:,2))
%% Finding out the frequency at each zero crossings

load('voltage.mat');
zero=circshift(voltage(:,2),1);
zero=zero.*voltage(:,2);
zero=zero<0;
zero1=circshift(zero,-1);
x2=voltage(zero,1);
%x2=circshift(x2,-1);
y2=voltage(zero,2);
%y2=circshift(y2,-1);
x1=voltage(zero1,1);
y1=voltage(zero1,2);
numer=x2.*y1-x1.*y2;
denom=y1-y2;
time=numer./denom;
time_start=circshift(time,2);
time_start(1:2)=0;
period=time-time_start;
freq_detected=[1./period(2:end),time(2:end,1)];
freq_detected(:,1)=round(smooth(freq_detected(:,1),20),2);
freq_detected(:,1)=round(smooth(freq_detected(:,1),20),2);
freq_steps=linspace(start_freq,end_freq,((end_freq-start_freq)/freq_resolution)+1);
freq_map=zeros(length(freq_steps),3);
freq_map(:,1)=freq_steps';
toc;
plot(freq_detected(:,2),freq_detected(:,1));
%% Mapping the frequency with the corresponding time start and end
for i=1:length(freq_map)
    temp_idx=freq_detected(:,1)==start_freq+((i-1)*freq_resolution);
    temp=freq_detected(temp_idx,2);
    if length(temp)>1
        first=min(temp);
        second=max(temp);
        temp_idx1=1;
        %deciding how much one frequency resolution time can be
        while(second-first>5)
            second=temp(length(temp)-temp_idx1);
            temp_idx1=temp_idx1+1;
        end
        freq_map(i,2)=first;
        freq_map(i,3)=second;
    end
end
plot(freq_map(:,1),freq_map(:,2),freq_map(:,1),freq_map(:,3));
%% Formulating the time vector for every frequency step
max_points=(round((freq_map(1,3)-freq_map(1,2)))+1)*sampling_rate;
%constructing force matrix
force_matrix=zeros(max_points,2*length(freq_map));
for j=1:length(freq_map)
    temp_idx=voltage(:,1)>freq_map(j,2)&voltage(:,1)<freq_map(j,3);
    temp_matrix=force(temp_idx,:);
    force_matrix(1:length(temp_matrix),(j*2)-1:j*2)=temp_matrix;
end
toc;
%constructing displacement matrix
disp_matrix=zeros(max_points,2*length(freq_map));
for j=1:length(freq_map)
    temp_idx=voltage(:,1)>freq_map(j,2)&voltage(:,1)<freq_map(j,3);
    temp_matrix=disp(temp_idx,:);
    disp_matrix(1:length(temp_matrix),(j*2)-1:j*2)=temp_matrix;
end
toc;
%save('E:\Thesis\aorta\test_cases\Rubber_Panel\2N_UP\disp_matrix.mat','disp_matrix')
%save('E:\Thesis\aorta\test_cases\Rubber_Panel\2N_UP\force_matrix.mat','force_matrix')

%%
results=zeros(length(freq_map),13);
for i=1:length(freq_map)
    X=disp_matrix(:,(i*2)-1:(i*2));
    Y=force_matrix(:,(i*2)-1:(i*2));
%     plot(X(:,2),Y(:,2));
%     pause(5);
%     close all;
    [area,max_Force,max_Disp,kk,coeff]=fitted_values(X,Y);
    results(i,1)=freq_map(i);
    results(i,2)=area;
    results(i,3)=area/(0.5*max_Force*max_Disp);
    results(i,4:7)=kk';
    results(i,8)=max_Force;
    results(i,9)=max_Disp;
    results(i,10:11)=coeff(2:3)';
    results(i,12:13)=coeff(6:7)';
end
k_fit=fit(results(:,1).^2,results(:,5),'poly1');
k_val=coeffvalues(k_fit);
figure;
subplot(2,2,1);
plot(freq_map(:,1),results(:,2));
xlabel('Frequency, Hz');
ylabel('Amplitude');
subplot(2,2,2);
plot(freq_map(:,1),results(:,2)./(pi*k_val(2)*results(:,9).^2));
xlabel('Frequency, Hz');
ylabel('K1, N/mm');
subplot(2,2,3);
plot(freq_map(:,1),results(:,5));
xlabel('Frequency, Hz');
ylabel('K2, N/mm^2');
subplot(2,2,4);
plot(freq_map(:,1),results(:,6));
xlabel('Frequency, Hz');
ylabel('K3, N/mm^3');
toc;
% %% 
set(0, 'DefaultFigureColor', 'White', ...
'DefaultFigurePaperType', 'a4letter', ...
'DefaultAxesColor', 'white', ...
'DefaultAxesFontUnits', 'points', ...
'DefaultAxesFontSize', 26, ...
'DefaultAxesFontAngle', 'normal', ...
'DefaultAxesFontName', 'Times New Roman', ...
'DefaultAxesGridLineStyle', ':', ...
'DefaultAxesInterruptible', 'on', ...
'DefaultAxesLayer', 'Bottom', ...
'DefaultAxesUnits', 'normalized', ...
'DefaultAxesXcolor', [0, 0, 0], ...
'DefaultAxesYcolor', [0, 0, 0], ...
'DefaultAxesZcolor', [0, 0, 0], ...
'DefaultAxesVisible', 'on', ...
'DefaultLineColor', 'Red', ...
'DefaultLineLineStyle', '-', ...
'DefaultLineLineWidth', 2, ...
'DefaultLineMarker', 'none', ...
'DefaultLineMarkerSize', 8, ...
'DefaultTextColor', [0, 0, 0], ...
'DefaultTextFontUnits', 'Points', ...
'DefaultTextFontSize', 30, ...
'DefaultTextFontName', 'Times', ...
'DefaultTextVerticalAlignment', 'middle', ...
'DefaultTextHorizontalAlignment', 'left');
%% 
% clear all;
% disp=load('disp.mat');
% force=load('force.mat');
% %removing the zeros and cleaning the data
% disp=disp.disp;
% force=force.force;
% X=disp;
% Y=[force(:,1),force(:,2)];
% temp_idx=X(:,1)~=0;
% X=X(temp_idx,:);
% temp_idx=Y(:,1)~=0;
% Y=Y(temp_idx,:);
% %since the force is having a negative sign
% Y=[Y(:,1),-Y(:,2)];
% 
% %plotting one cycle
% temp=circshift(X(:,2),1);
% temp(1)=0;
% temp_idx=(X(:,2).*temp)<0;
% temp_idx1=find(temp_idx);
% X=X(temp_idx1(1):temp_idx1(3),:);
% Y=Y(temp_idx1(1):temp_idx1(3),:);
% [val1,minidx]=min(X(:,2));
% [val2,maxidx]=max(X(:,2));
% val3=min(Y(:,2));
% %possible confusion in index of minimum and maximum values
% x1=X(maxidx:minidx,2);
% y1=Y(maxidx:minidx,2);
% x2=X(1:maxidx,2);
% x3=X(minidx:end,2);
% x2=[x2;x3];
% y2=Y(1:maxidx,2);
% y3=Y(minidx:end,2);
% y2=[y2;y3];
% plot(x1,y1,x2,y2);
% x3=x1-val1;
% y3=y1-val3;
% x4=x2-val1;
% y4=y2-val3;
% [values,order]=sort(x4);
% x4=x4(order,:);
% y4=y4(order,:);
% [values,order]=sort(x3);
% x3=x3(order,:);
% y3=y3(order,:);
% area=abs(trapz(x3,y3)-trapz(x4,y4));
% 
% %% 
% xdMat=zeros(4656,2);
% ydMat=zeros(4656,2);
% energy=zeros(length(freq_map),7);
% for k=26:50
%     idx=force_matrix(:,k*2)~=0;
%     force_temp=force_matrix(idx,k*2-1:k*2);
%     disp_temp=disp_matrix(idx,k*2-1:k*2);
%     no_of_points=round(length(force_temp)-(1/freq_map(k,1))*4*sampling_rate);
%     xdMat(1:length(disp_temp(no_of_points:end,2)),k-1)=disp_temp(no_of_points:end,2);
%     ydMat(1:length(disp_temp(no_of_points:end,2)),k-1)=force_temp(no_of_points:end,2);
%     XY=[disp_temp(no_of_points:end,2),force_temp(no_of_points:end,2)];
%     figure;
%     plot(XY(:,1),XY(:,2),'+');
%     hold on;
%     A=EllipseDirectFit(XY);
%     energy(k,1:6)=A';
%     %Convert the A to str 
%     a = num2str(A(1)); 
%     b = num2str(A(2)); 
%     c = num2str(A(3)); 
%     d = num2str(A(4)); 
%     e = num2str(A(5)); 
%     f = num2str(A(6));
%     aa = str2double(a);
%     bb = str2double(b)/2; 
%     cc = str2double(c); 
%     dd = str2double(d)/2; 
%     ff = str2double(e)/2; 
%     gg = str2double(f);
%      eqt= ['(',a, ')*x^2 + (',b,')*x*y + (',c,')*y^2 + (',d,')*x+ (',e,')*y + (',f,')']; 
%     xmin=-0.06; 
%     xmax=0.06; 
%     ezplot(eqt,[xmin,xmax]) 
%     num=2*((aa*ff^2+cc*dd^2+gg*bb^2-(2*bb*dd*ff)-(aa*cc*gg)));
%     denA=(bb^2-aa*cc)*(sqrt((aa-cc)^2+4*bb^2)-(aa+cc));
%     denB=(bb^2-aa*cc)*(-sqrt((aa-cc)^2+4*bb^2)-(aa+cc));
%     adash=sqrt(num/denA);
%     bdash=sqrt(num/denB);
%     area=adash*bdash*pi;
%     energy(k,7)=area;
%     hold off;
% end
% 
% %% 
% 
% zMat = freq_map(2:44,1); %// For plot3
% 
% plot3(xdMat,zMat, ydMat, '+'); %// Make all traces blue
% grid;
% xlabel('x'); ylabel('y'); zlabel('z'); %// Adjust viewing angle so you can clearly see data
% %% 
% %% 
% plot(disp_matrix(30000:33683,20),force_matrix(30000:33683,20))
% hold on
% XY=[disp_matrix(30000:33683,20),force_matrix(30000:33683,20)];
% A=EllipseDirectFit(XY);
% %Convert the A to str 
% a = num2str(A(1)); 
% b = num2str(A(2)); 
% c = num2str(A(3)); 
% d = num2str(A(4)); 
% e = num2str(A(5)); 
% f = num2str(A(6));
% aa = str2double(a);
% bb = str2double(b)/2; 
% cc = str2double(c); 
% dd = str2double(d)/2; 
% ff = str2double(e)/2; 
% gg = str2double(f);
% %Equation 
% eqt= ['(',a, ')*x^2 + (',b,')*x*y + (',c,')*y^2 + (',d,')*x+ (',e,')*y + (',f,')']; 
% xmin=-0.06; 
% xmax=0.06; 
% ezplot(eqt,[xmin,xmax]) 
% num=2*((aa*ff^2+cc*dd^2+gg*bb^2-(2*bb*dd*ff)-(aa*cc*gg)));
% denA=(bb^2-aa*cc)*(sqrt((aa-cc)^2+4*bb^2)-(aa+cc));
% denB=(bb^2-aa*cc)*(-sqrt((aa-cc)^2+4*bb^2)-(aa+cc));
% adash=sqrt(num/denA);
% bdash=sqrt(num/denB);
% area=adash*bdash*pi
% hold off