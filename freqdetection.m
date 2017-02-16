clear all;
close all;
tic;
%Constants
freq_resolution=1;
start_freq=10;
end_freq=100;
sampling_rate=12800;
start_time=4;
end_time=120;
no_of_periods=4;

voltage =  importdata('E:\Thesis\aorta\100mmHg\100mmhg_voltage.txt','\t',31);
force =  importdata('E:\Thesis\aorta\100mmHg\100mmhg_force.txt','\t',31);
disp =  importdata('E:\Thesis\aorta\100mmHg\100mmhg_disp.txt','\t',31);
voltage = voltage.data;

start_time_idx=voltage(:,1)>start_time & voltage(:,1)<end_time;
voltage=voltage(start_time_idx,:);
force=force.data(start_time_idx,:);
disp=disp.data(start_time_idx,:);
save('voltage.mat','voltage');
save('force.mat','force');
save('disp.mat','disp');
toc;
% plot(voltage(:,1),voltage(:,2))
% Finding out the frequency at each zero crossings

load('voltage.mat');
zero=circshift(voltage(:,2),1);
zero=zero.*voltage(:,2);
zero=zero<0;
time=voltage(zero,1);
time_start=circshift(time,1);
time_start(1)=0;
period=time-time_start;
freq_detected=[round(1./period(2:end)),time(2:end,1)];
freq_steps=linspace(start_freq,end_freq,((end_freq-start_freq)*freq_resolution)+1);
freq_map=zeros(length(freq_steps),3);
freq_map(:,1)=freq_steps';
toc;
% Mapping the frequency with the corresponding time start and end

freq_count=1;
for i=1:length(freq_detected)
   if i==1
       freq_map(freq_count,1:2)=freq_detected(1,:);
   else
       if freq_detected(i,1)==freq_detected(i-1,1)+freq_resolution
           freq_map(freq_count,3)=freq_detected(i-1,2);
           freq_count=freq_count+1;
           freq_map(freq_count,1:2)=freq_detected(i,:);
       end
   end
end
% Formulating the time vector for every frequency step
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
% 
xdMat=zeros(4656,2);
ydMat=zeros(4656,2);
energy=zeros(length(freq_map),7);
for k=26:50
    idx=force_matrix(:,k*2)~=0;
    force_temp=force_matrix(idx,k*2-1:k*2);
    disp_temp=disp_matrix(idx,k*2-1:k*2);
    no_of_points=round(length(force_temp)-(1/freq_map(k,1))*4*sampling_rate);
    xdMat(1:length(disp_temp(no_of_points:end,2)),k-1)=disp_temp(no_of_points:end,2);
    ydMat(1:length(disp_temp(no_of_points:end,2)),k-1)=force_temp(no_of_points:end,2);
    XY=[disp_temp(no_of_points:end,2),force_temp(no_of_points:end,2)];
    figure;
    plot(XY(:,1),XY(:,2),'+');
    hold on;
    A=EllipseDirectFit(XY);
    energy(k,1:6)=A';
    %Convert the A to str 
    a = num2str(A(1)); 
    b = num2str(A(2)); 
    c = num2str(A(3)); 
    d = num2str(A(4)); 
    e = num2str(A(5)); 
    f = num2str(A(6));
    aa = str2double(a);
    bb = str2double(b)/2; 
    cc = str2double(c); 
    dd = str2double(d)/2; 
    ff = str2double(e)/2; 
    gg = str2double(f);
     eqt= ['(',a, ')*x^2 + (',b,')*x*y + (',c,')*y^2 + (',d,')*x+ (',e,')*y + (',f,')']; 
    xmin=-0.06; 
    xmax=0.06; 
    ezplot(eqt,[xmin,xmax]) 
    num=2*((aa*ff^2+cc*dd^2+gg*bb^2-(2*bb*dd*ff)-(aa*cc*gg)));
    denA=(bb^2-aa*cc)*(sqrt((aa-cc)^2+4*bb^2)-(aa+cc));
    denB=(bb^2-aa*cc)*(-sqrt((aa-cc)^2+4*bb^2)-(aa+cc));
    adash=sqrt(num/denA);
    bdash=sqrt(num/denB);
    area=adash*bdash*pi;
    energy(k,7)=area;
    hold off;
end

% 

zMat = freq_map(2:44,1); %// For plot3

plot3(xdMat,zMat, ydMat, '+'); %// Make all traces blue
grid;
xlabel('x'); ylabel('y'); zlabel('z'); %// Adjust viewing angle so you can clearly see data%% 
% 
plot(disp_matrix(30000:33683,20),force_matrix(30000:33683,20))
hold on
XY=[disp_matrix(30000:33683,20),force_matrix(30000:33683,20)];
A=EllipseDirectFit(XY);
%Convert the A to str 
a = num2str(A(1)); 
b = num2str(A(2)); 
c = num2str(A(3)); 
d = num2str(A(4)); 
e = num2str(A(5)); 
f = num2str(A(6));
aa = str2double(a);
bb = str2double(b)/2; 
cc = str2double(c); 
dd = str2double(d)/2; 
ff = str2double(e)/2; 
gg = str2double(f);
%Equation 
eqt= ['(',a, ')*x^2 + (',b,')*x*y + (',c,')*y^2 + (',d,')*x+ (',e,')*y + (',f,')']; 
xmin=-0.06; 
xmax=0.06; 
ezplot(eqt,[xmin,xmax]) 
num=2*((aa*ff^2+cc*dd^2+gg*bb^2-(2*bb*dd*ff)-(aa*cc*gg)));
denA=(bb^2-aa*cc)*(sqrt((aa-cc)^2+4*bb^2)-(aa+cc));
denB=(bb^2-aa*cc)*(-sqrt((aa-cc)^2+4*bb^2)-(aa+cc));
adash=sqrt(num/denA);
bdash=sqrt(num/denB);
area=adash*bdash*pi
hold off