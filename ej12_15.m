clear all ; close all; clc
T=2*pi; 

t = -T/2:1e-4:T/2;

k=100;
sf=0;

f= -(t/pi).^2.*(t>=-T/2 & t <0) + ((t/pi).^2.*(t>=0 & t <T/2));
plot(t, f, 'r'), hold 'on'

for n= 1:k
    bn = 2*((2-pi^2*n^2)*cos(n*pi)-2)/(pi^3*n^3);
    sf = sf + bn*sin(n*t);

endfor

plot(t,sf)
  
