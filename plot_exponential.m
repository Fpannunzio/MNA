clear all ; close all; clc
A=1; T=1; 

t = 0:1e-4:T;

n=150;
sf=0;

f= A.*(t>0 & t <=T/2) - A.*(t>T/2 & t<1);
plot(t, f, 'r'), hold 'on'

for k=-n:n
  if k == 0
    continue
  endif
  cn = (A/(i*pi*k))*(1-e.^(-i*pi*k));
  sf = sf + cn*(e.^(i*2*pi*k*t/T));
end

plot(t,sf)