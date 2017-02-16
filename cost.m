function J = cost(theta1)
global B A
m = length(B); % number of training examples
predictions=A*theta1';
weight=ones(length(B),1);
sqrerror=((predictions-B).*weight).^2;
J=(1/(2*m))*(sum(sqrerror));
end