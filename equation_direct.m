function dy = equation_direct(t,q)
global f m omega zeta k k2 k3
dy =zeros(2,1);
dy(1)=q(2);
dy(2)=(f/m).*(cos(t.*omega))-(k/m).*q(1)-2*zeta*sqrt(k/m).*q(2)-(k2/m).*q(1).^2-(k3/m).*q(1).^3;
end