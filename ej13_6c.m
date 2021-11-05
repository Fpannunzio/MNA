clear all ; close all; clc
T=5/127; 

t = 1;

fk = []

N=127;

for k= 0:N
  fk=0;
  for n= 0:N
    cn = e^-(n*T);
    fk(k) = fk(k) + cn*e.^(i*2*pi*k*n/(N+1));
  endfor
  fk(k)
endfor

