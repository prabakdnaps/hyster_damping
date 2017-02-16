clear all;
close all;
tic;
fold_name='E:\Thesis\aorta\test_cases\Rubber_Panel\0.5N_UP\';
volt_filename='0.5N_RP_up_force.txt';
force_filename='0.5N_RP_up_force.txt';
disp_filename='0.5N_RP_up_disp.txt';

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

% freq_resolution=0.05;
% start_freq=82;
% end_freq=110;
% sampling_rate=6400;
% start_time=4;
% end_time=970;
% no_of_periods=4;
% force =  importdata('E:\Thesis\aorta\test_cases\Rubber_Panel\2N_UP\2N_UP_force.txt','\t',34);
% disp =  importdata('E:\Thesis\aorta\test_cases\Rubber_Panel\2N_UP\2N_UP_disp.txt','\t',34);

freq_resolution=0.05;
start_freq=98;
end_freq=108;
sampling_rate=6400;
start_time=4;
end_time=135;
no_of_periods=4;
force =  importdata(strcat(fold_name,force_filename),'\t',34);
disp =  importdata(strcat(fold_name,disp_filename),'\t',34);


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
figure;
subplot(2,2,1);
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
subplot(2,2,2);
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
subplot(2,2,3);
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
save(strcat(fold_name,'disp_matrix.mat'),'disp_matrix')
save(strcat(fold_name,'force_matrix.mat'),'force_matrix')

%%

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
