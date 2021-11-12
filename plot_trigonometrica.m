clear all ; close all; clc
A=1; T=1; 

t = 0:1e-4:T;

n=150;
sf=0;

f=(-1/2 + 2*t).*(t>0 & t <=1/2) +(3/2 -2*t).*(t>1/2 & t<1);
plot(t, f, 'r'), hold 'on'

for k=1:n
  bn= 0;
  an= 2*(cos(pi*k) -1)/((pi**2)*(k**2));
  sf = sf + an*cos(2*pi*k*t) + bn*sin(2*pi*k*t);
end

plot(t,sf)