clear all; clc


h=0.1;
tf=1000;
y0=0.35503;

x = -2 - [0:h:tf]*h^2;
n=size(x)(2)-1

M = diag(ones(n,1),1) + diag(ones(n,1),-1) + diag(x);

V = [ -y0 zeros(1, n)];

P = inv(M)*V';

P(size(P)(1)) 
airy(0,tf)