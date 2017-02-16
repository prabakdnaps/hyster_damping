global f m omega zeta k k2 k3
tic;
f=2.0;
test=linspace(80,110,20);
num_results=zeros(10,6);
for i=1:length(test)
    omega=test(i)*2*pi;
    zeta =0.05;
    k2=10.0;
    k3=-6000;
    k=300;
    m=7e-4;
    t=linspace(0,100,5000);
    y0=[0,0];
    [T,Y]=ode45(@equation_direct,t,y0);
    force=f.*cos(t.*omega);
    test_case1=[T,Y(:,1),force'];
    idx=test_case1(:,1)>50;
    test_case1=test_case1(idx,:);
    X=test_case1(:,1:2);
    Y=[test_case1(:,1),test_case1(:,3)];
    num_results(i,1)=test(i);
    num_results(i,2)=max(X(:,2));
    toc;
end