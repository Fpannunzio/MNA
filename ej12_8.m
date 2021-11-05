clear all ; close all; clc
A=1; T=1; 

t = 0:1e-4:T;

k=100;
sf=0;

f= 1.*(t>0 & t <T/2) + -1.*(t>T/2 & t<T);
plot(t, f, 'r'), hold 'on'

for n= 1:k
  n, cn = (i*e^(-2*i*pi*n)*(-1 + e^(i* pi* n)))/(2*  n) -(i - i *e^(-i* pi *n))/(2* pi *n)
  sf = sf + cn*e.^(i*2*pi*n*t/T);
end

#for n= -k:-1
  #n, cn = (i*e^(-2*i*pi*n)*(-1 + e^(i* pi* n)))/(2*  n) -(i - i *e^(-i* pi *n))/(2* pi *n)
  #sf = sf + cn*e.^(-i*2*pi*n*t/T);
#end

plot(t,sf)