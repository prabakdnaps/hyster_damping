A=zeros(2*length(freq_map),1);
B=ones(2*length(freq_map),1);
wn=103.25;
zeta=0.035;
for i=1:length(freq_map)
A(2*i-1,1)=(-(results(i,1)/wn)^2*results(i,10)+2*(results(i,1)/wn)*zeta*results(i,11)+results(i,10))/results(i,12);
A(2*i,1)=(-(results(i,1)/wn)^2*results(i,11)-2*(results(i,1)/wn)*zeta*results(i,10)+results(i,11))/results(i,13);
end
options = optimoptions(@lsqlin,'Display','none');
k_est = lsqlin(A,B,[],[],[],[],[0;-Inf;-Inf],[],[],options);
%% 

k=75.5*5;
m=0.000891;
c=2*m*0.035*103.25*2*pi;
omega_pre=2*pi*results(:,1);
X_est=3.2./((k-(m.*omega_pre.^2))+1i*c.*omega_pre);
plot(omega_pre,abs(X_est));
hold on
plot(results(:,1).*2*pi,results(:,9))