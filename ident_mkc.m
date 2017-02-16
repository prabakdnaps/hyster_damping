function [kk] = ident_mkc(X,Y)

coeff=zeros(1,8);
foco1=fit(X(:,1),X(:,2),'fourier1');
coeff(1:4)=coeffvalues(foco1);
foco2=fit(Y(:,1),Y(:,2),'fourier1');
coeff(5:8)=coeffvalues(foco2);
mm=coeff(4).^2*(-coeff(2).*cos(coeff(4).*X(:,1))-coeff(3).*sin(coeff(4).*X(:,1)));
cm=coeff(4)*(-coeff(2).*sin(coeff(4).*X(:,1))+coeff(3).*cos(coeff(4).*X(:,1)));
km=coeff(1)+coeff(2).*cos(coeff(4).*X(:,1))+coeff(3).*sin(coeff(4).*X(:,1));
fm=coeff(6).*cos(coeff(8).*X(:,1))+coeff(7).*sin(coeff(8).*X(:,1));
k2=coeff(2).^2.*(cos(coeff(4).*X(:,1))).^2+coeff(3).^2.*(sin(coeff(4).*X(:,1))).^2+2*coeff(2)*coeff(3).*cos(coeff(4).*X(:,1)).*sin(coeff(4).*X(:,1));
A=[mm,cm,km];
b=fm;
omega_n=1;
zeta=0.05;
k1m=mm./(omega_n).^2;
k2m=(2*zeta/omega_n).*cm;
k3m=km;
Ak=[km.^2,km.^3];
kk=Ak\(-(k1m+k2m+k3m)+fm);
end
