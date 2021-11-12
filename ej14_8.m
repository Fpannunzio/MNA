clear all; clc


h=0.1;
tf=1;
y0=1;
a = 2;

n=1/h

x = ones(n+1,1) * (- 2 - a^2*h^2);

M = diag(ones(n,1),1) + diag(ones(n,1),-1) + diag(x);

M(n+1,n) = 2

V = [ -y0 zeros(1, n-1) -2*h];

P = inv(M)*V';

P(size(P)(1)) 
cosh(0)/cosh(a)