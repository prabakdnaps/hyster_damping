theta = [1;1;1];% initialize fitting parameters
% Some gradient descent settings
iterations = 200000;
alpha = 0.01;
% run gradient descent
[theta,J] = gradientDescentMulti(A,B, theta, alpha, iterations);
plot(J);
% print theta to screen
fprintf('Theta found by gradient descent: ');
theta