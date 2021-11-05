clear all ; close all; clc
T=1; 

t = 1;

k=100;
sf=0;

f= e.^(-5*t).*(t>=0 & t <T);
plot(t, f, 'r'), hold 'on'

for n= -k:k
  cn = (1-e^-(5+2*pi*n*i))/(5 + 2*pi*i*n);
  sf = sf + cn*e.^(i*2*pi*n*t);
endfor

sf

plot(t,sf)
  
