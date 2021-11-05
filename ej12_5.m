clear all ; close all; clc
A=1; T=1; 

t = 0:1e-2:T;

n=50;
sf=0;

f=(-1/2 + 2*t).*(t>0 & t <=1/2) +(3/2 -2*t).*(t>1/2 & t<1);
plot(t, f, 'r'), hold 'on'

for k=1:n
  bn= 2*((cos((3*pi*k)/2)*(pi*k*cos((pi*k)/2)))/(2*pi^2*k^2) - (pi* k + pi* k *cos(pi* k))/(4*pi^2 *k^2));
  an= 2*((2*cos(pi*k)-2)/(4*k^2*pi^2) - (sin(3*pi*k/2)*(pi*k*cos(pi*k/2)-2*sin(pi*k/2)))/(2*pi^2*k^2));
  sf = sf + an*cos(2*pi*k*t) + bn*sin(2*pi*k*t);
end

plot(t,sf)