clear all;
N1=load('E:\Thesis\aorta\test_cases\Rubber_Panel\0.1N\results.mat');
N5=load('E:\Thesis\aorta\test_cases\Rubber_Panel\0.5N_UP\results.mat');
area1=smooth(N1.results(:,2),21);
store1=smooth(N1.results(:,9).^2,21);
area2=smooth(N5.results(:,2),21);
store2=smooth(N5.results(:,9).^2,21);
loss1=smooth(area1./(pi*306*store1),21);
loss2=smooth(area2./(pi*306*store2),21);
plot(N1.results(:,1),loss1,'o');
hold on
plot(N5.results(:,1),loss2,'+');
zeta=zeros(length(loss2),1);
zeta=zeta+2*0.0388;
plot(N5.results(:,1),zeta)
legend('0.1N','0.5N','2*\zeta');
xlabel('Frequency, Hz');
ylabel('Loss Factor, \\');
figure;
plot(N1.results(:,1),area1,'o');
hold on
plot(N1.results(:,1),store1*pi*306,'+');
legend('Area under curve','Storage energy');
xlabel('Frequency, Hz');
ylabel('Energy,N.mm');
figure;
plot(N5.results(:,1),area2,'o');
hold on
plot(N5.results(:,1),store2*pi*306,'+');
legend('Area under curve','Storage energy');
xlabel('Frequency, Hz');
ylabel('Energy,N.mm');