clear all ; close all; clc
A=1; T=pi; 

t = -pi:1e-4:pi;

n=150;
sf=0;

f=(-1/2 + 2*t).*(t>0 & t <=1/2) +(3/2 -2*t).*(t>1/2 & t<1);
f=t;
plot(t, f, 'r'), hold 'on'

for k=1:n
  bn= 2*(-1)^(k+1)/k;
  an= 0;
  sf = sf + an*cos(k*t) + bn*sin(k*t);
end

plot(t,sf)