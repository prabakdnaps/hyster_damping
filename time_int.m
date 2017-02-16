%% test case1
global f ratio zeta k2 k3
tic;
f=0.2;
test=linspace(0.6,1.4,40);
results=zeros(10,6);
for i=1:length(test)
    ratio=test(i);
    zeta =0.05;
    k2=0;
    k3=0;
    t=linspace(0,500,50000);
    y0=[0,0];
    [T,Y]=ode45(@equation,t,y0);
    force=f.*cos(t.*ratio);
    test_case1=[T,Y(:,1),force'];
    idx=test_case1(:,1)>400;
    test_case1=test_case1(idx,:);
    X=test_case1(:,1:2);
    Y=[test_case1(:,1),test_case1(:,3)];
    [area,min_Disp,max_Disp,KK,kk]=fitted_values(X,Y);
%     results(i,1)=ratio;
%     results(i,2)=area;
%     results(i,3)=area/(0.5*max_Force*max_Disp);
%     results(i,4:7)=kk';
%     results(i,8)=max_Force;
%     results(i,9)=max_Disp;
%     results(i,10:11)=coeff(2:3)';
%     results(i,12:13)=coeff(6:7)';
   % kk=ident_mkc(X,Y);
    results(i,1)=ratio;
    results(i,2)=area;
    results(i,4)=max_Disp;
    results(i,5)=min_Disp;
    results(i,6:9)=KK';
    results(i,10)=kk(1);
    results(i,11)=kk(2);
    results(i,12)=kk(3);
end
toc;
k1=1;
k2=mean(results(:,11));
k3=mean(results(:,12));
for i=1:length(results)
results(i,3)=pi*(((k1/2*results(i,4).^2)+(k2/3*results(i,4).^3)+...
(k3/4*results(i,4).^4))+((k1*results(i,5).^2)-(k2/3*results(i,5).^3)...
+(k3/4*results(i,5).^4)));
end
figure;
subplot(2,2,1);
plot(results(:,1),results(:,2));
xlabel('Frequency, Hz');
ylabel('Amplitude');
subplot(2,2,2);
plot(results(:,1),results(:,2)./results(:,3));
xlabel('Frequency, Hz');
ylabel('K1, N/mm');
subplot(2,2,3);
plot(results(:,1),results(:,4));
xlabel('Frequency, Hz');
ylabel('K2, N/mm^2');
subplot(2,2,4);
plot(results(:,1),results(:,5));
xlabel('Frequency, Hz');
ylabel('K3, N/mm^3');
%% 
global f ratio zeta k2 k3
f=4;
ratio=0.89;
zeta =0.05;
k2=0;
k3=0.2;
t=linspace(0,500,50000);
y0=[0,0];
[T,Y]=ode45(@equation,t,y0);
force=f.*cos(t.*ratio);
test_case2=[T,Y(:,1),force'];
idx=test_case2(:,1)>400;
test_case2=test_case2(idx,:);
save('test_case2.mat','test_case2');
figure;
subplot(1,2,1);
plot(test_case2(:,1),test_case2(:,2),test_case2(:,1),test_case2(:,3));
xlabel('Time, s');
ylabel('Amplitude, N & mm');
subplot(1,2,2);
plot(test_case2(:,2),test_case2(:,3))
xlabel('Displacement, mm');
ylabel('Force, N');

%% 
global f ratio zeta k2 k3
f=0.2;
ratio=0.9;
zeta =0.05;
k2=0.2;
k3=0;
t=linspace(0,500,50000);
y0=[0,0];
[T,Y]=ode45(@equation,t,y0);
force=f.*cos(t.*ratio);
test_case3=[T,Y(:,1),force'];
idx=test_case3(:,1)>400;
test_case3=test_case3(idx,:);
save('test_case3.mat','test_case3');
figure;
subplot(1,2,1);
plot(test_case3(:,1),test_case3(:,2),test_case3(:,1),test_case3(:,3));
xlabel('Time, s');
ylabel('Amplitude, N & mm');
subplot(1,2,2);
plot(test_case3(:,2),test_case3(:,3))
xlabel('Displacement, mm');
ylabel('Force, N');

%% 

global f ratio zeta k2 k3
f=0.2;
ratio=0.99;
zeta =0.02;
k2=0;
k3=0;
t=linspace(0,500,50000);
y0=[0,0];
[T,Y]=ode45(@equation,t,y0);
force=f.*cos(t.*ratio);
test_case3=[T,Y(:,1),force'];
idx=test_case3(:,1)>400;
test_case3=test_case3(idx,:);
save('test_case3.mat','test_case3');
% figure;
subplot(1,2,1);
plot(test_case3(:,1),test_case3(:,2),test_case3(:,1),test_case3(:,3));
xlabel('Time, s');
ylabel('Amplitude, N & mm');
subplot(1,2,2);
plot(test_case3(:,2),test_case3(:,3))
xlabel('Displacement, mm');
ylabel('Force, N');


