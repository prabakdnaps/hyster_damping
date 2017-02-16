hz=linspace(75,125,500);
omega=hz.*2*pi;
f=0.5;
k=300;
m=k/(103.5*2*pi)^2;
zeta=0.038;
c=2*m.*sqrt(k/m)*zeta;
x=f./(sqrt((k-m.*omega.^2).^2+(c.*omega).^2));
plot(hz,x)
% y=0.5*ones(length(x),1)-(1i*omega.*x)';
% omega_ratio=omega./(103.25*2*pi);
% X=(x-omega_ratio.^2.*x)';
% act=[k;m;c];
% X(:,1)=X(:,1)*act(1);
% X(:,2)=X(:,2)*act(2);
% theta=rand(2,1);
% alpha=0.001;
% lambda=0;
% num_iters=10000;
% % X_mean=mean(X);
% % X_std=std(X);
% % for i=1:3
% %     X(:,i)=(X(:,i))./X_std(i);
% % end
% [Theta,J]=gradientDescentMulti(X, y, theta, alpha,lambda, num_iters);
% plot(abs(J))
%% 
freq_hz=results(:,1);
freq_ratio=results(:,1)./103.25;
freq_ratio2=freq_ratio.^2;
X1=-freq_ratio2.*results(:,10)+2*0.035.*freq_ratio.*results(:,11)+results(:,10);
Y1=results(:,12);
X2=-freq_ratio2.*results(:,11)-2*0.035.*freq_ratio.*results(:,10)+results(:,11);
Y2=results(:,13);
XX=zeros(2*length(freq_hz),1);
YY=zeros(2*length(freq_hz),1);
for i=1:length(freq_hz)
    XX(2*i-1)=X1(i);
    XX(2*i)=X2(i);
    YY(2*i-1)=Y1(i);
    YY(2*i)=Y2(i);
end
XX\YY
theta=0;
alpha=0.001;
lambda=0;
num_iters=10000;
XX=XX./0.0016;
YY=YY./0.5;
[Theta,J]=gradientDescentMulti(XX, YY, theta, alpha,lambda, num_iters);
plot(J)