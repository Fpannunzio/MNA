clear all ; close all; clc
T=1; 

t = 0:1e-3:T;

k=100;
sf=0;

f= e.^(-5*t).*(t>=0 & t <T);
#plot(t, f, 'r'), hold 'on'

for n= -k:k
  
  if n != 0
    cn = (1-e^-(pi*n*i))/(2*pi*i*n);
    sf = sf + cn*e.^(i*2*pi*n*t);
  else
    sf = sf + 1/2;
  endif
endfor

plot(t,sf)
  
