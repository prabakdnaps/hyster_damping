function dy = equation(t,q)
global f ratio zeta k2 k3
dy =zeros(2,1);
dy(1)=q(2);
dy(2)=f.*(cos(t.*ratio))-q(1)-2*zeta*q(2)-k2*q(1).^2-k3*q(1).^3;
end