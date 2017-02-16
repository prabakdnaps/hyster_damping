%folder name in which the analysis to be done
fold_name='E:\Thesis\aorta\test_cases\Rubber_Panel\0.5N_UP\';
load(strcat(fold_name,'disp_matrix.mat'))
load(strcat(fold_name,'disp_matrix.mat'))

tic;
results=zeros(length(freq_map),13);
for i=1:length(freq_map)
    X=disp_matrix(:,(i*2)-1:(i*2));
    Y=force_matrix(:,(i*2)-1:(i*2));
    [area,min_Disp,max_Disp,KK,kk]=fitted_values(X,Y);
    results(i,1)=freq_map(i);
    results(i,2)=area;
    results(i,5)=max_Disp;
    results(i,6)=min_Disp;
    results(i,7:10)=KK';
    results(i,11)=kk(1);
    results(i,12)=kk(2);
    results(i,13)=kk(3);
end

%identification of k1 from the slopes of the force vs displacement curves
% k_fit=fit(results(:,1).^2,results(:,5),'poly1');
% k_val=coeffvalues(k_fit);
% results(:,3)=pi*k_val(2)*results(:,9).^2;

%given k1, k2 and k3 values
k1=300;
k2=-0.0625;
k3=-3000;

%storage energy calculation
for i=1:length(freq_map)
results(i,3)=pi*(((k1/2*results(i,8).^2)+(k2/3*results(i,8).^3)+...
(k3/4*results(i,8).^4))+((k1/2*results(i,8).^2)-(k2/3*results(i,8).^3)...
+(k3/4*results(i,8).^4)));
end

%plotting the results

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
plot(freq_map(:,1),results(:,4));
xlabel('Frequency, Hz');
ylabel('Loss Factor');
subplot(2,2,4);
plot(freq_map(:,1),results(:,6));
xlabel('Frequency, Hz');
ylabel('K3, N/mm^3');
toc;