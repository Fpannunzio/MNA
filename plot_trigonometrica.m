clear all ; close all; clc
A=1; T=pi; 

t = -pi:1e-4:pi;

n=150;
sf=-pi^2/3;

f=(-1/2 + 2*t).*(t>0 & t <=1/2) +(3/2 -2*t).*(t>1/2 & t<1);
f=t.*(pi-t);
plot(t, f, 'r'), hold 'on'

for k=1:n
  bn= -2*pi*(-1)^(k)/k;
  an= -4*(-1)^(k)/k^2;
  sf = sf + an*cos(k*t) + bn*sin(k*t);
end

plot(t,sf)