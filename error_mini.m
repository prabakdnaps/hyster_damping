% x = [1000, 2000, 4000];
% y = [200000, 250000, 300000];
% X=[24.1103,7658.5;100.543,31938.2];
% y=[13424;66363.2];
% given a theta_0 and theta_1, this function calculates
% their cost. We don't need this function, strictly speaking...
% but it is nice to print out the costs as gradient descent iterates.
% We should see the cost go down every time the values of theta get updated.
options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton');
theta = (rand())*[1;1];
opt_theta = fminunc(@cost, theta,options);
error=A*opt_theta';
disp(error);
disp(opt_theta');