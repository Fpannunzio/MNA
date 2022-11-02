clear all ; close all; clc
A=1; T=2*pi; 

t = -pi:1e-2:pi;

n=15;
sf=0;

f= t.^2;
plot(t, f, 'r'), hold 'on'

for k=0:n
  if k == 0
    continue
  endif
  cn = 2*(-1)^(k)/k^2;
  sf = sf + cn*(e.^(i*k*t));
end

plot(t,sf)