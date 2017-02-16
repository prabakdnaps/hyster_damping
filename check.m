% load('disp_full.mat');
% load('force_full.mat');
tic;
results=zeros(length(freq_map),13);
for i=1:length(freq_map)
    X=disp_matrix(:,(i*2)-1:(i*2));
    Y=force_matrix(:,(i*2)-1:(i*2));
%     plot(X(:,2),Y(:,2));
%     pause(5);
%     close all;
    [area,min_Disp,max_Disp,KK,kk]=fitted_values(X,Y);
    results(i,1)=freq_map(i);
    results(i,2)=area;
    results(i,4)=max_Disp;
    results(i,5)=min_Disp;
    results(i,6:9)=KK';
    results(i,10)=kk(1);
    results(i,11)=kk(2);
    results(i,12)=kk(3);
end
k_fit=fit(results(:,1).^2,results(:,5),'poly1');
k_val=coeffvalues(k_fit);
k1=300;
k2=-0.0625;
k3=-3000;

for i=1:length(freq_map)
results(i,3)=pi*(((k1/2*results(i,8).^2)+(k2/3*results(i,8).^3)+...
(k3/4*results(i,8).^4))+((k1/2*results(i,8).^2)-(k2/3*results(i,8).^3)...
+(k3/4*results(i,8).^4)));
end
% results(:,3)=pi*k_val(2)*results(:,9).^2;
figure;
subplot(2,2,1);
plot(freq_map(:,1),results(:,2));
xlabel('Frequency, Hz');
ylabel('Area Under Curve');
subplot(2,2,2);
plot(freq_map(:,1),results(:,3));
xlabel('Frequency, Hz');
ylabel('Storage Energy');
subplot(2,2,3);
plot(freq_map(:,1),results(:,2)./results(:,3));
xlabel('Frequency, Hz');
ylabel('Loss Factor');
subplot(2,2,4);
plot(freq_map(:,1),results(:,6));
xlabel('Frequency, Hz');
ylabel('K3, N/mm^3');
toc;
%% test cases
load('test_case1.mat')
X=test_case1(:,1:2);
Y=[test_case1(:,1),test_case1(:,3)];
% nonl =  importdata('E:\Thesis\aorta\soft.txt','\t');
% plot(nonl(:,1),nonl(:,2))
% X=zeros(length(nonl),2);
% X(:,2)=nonl(:,1);
% Y=zeros(length(nonl),2);
% Y(:,2)=nonl(:,2);
[area,KK]=fitted_values(X,Y);