clear all; clc



tf=1000;
n=4;

h=1/(n+1)

x = ones(n-1,1)*(1-h^2);

M = diag(ones(n,1)*-2) + diag(ones(n-1,1),-1) + diag(x,1);

M(1,2) = (2-h^2);

V = [ ones(1, n-1)*h^2 (2*h^2 -1)];

P = inv(M)*V';

x=1;

f = 2*(e^x+e^(-x))/(e+e^(-1))-1;

ex = @(x) 2*(exp(x) + exp(-x))/(e + e^-1) - 1;


P(size(P)(1))

0.2439/ERR

