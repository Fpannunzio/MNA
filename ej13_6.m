clear all; close all; clc

T1 = 5/127;
F1 = 1/T1;
tfin = 5;

n = 0:T1:tfin;
L = length(n);
NP = 1;

xn = cos(n)

NFFT = 2^nextpow2(NP*L);
xf1 = fft(xn, NFFT)/L

k = 0:1:127;

f1 = (F1/2)*linspace(0, 1, NFFT/2);
stem(f1,2*abs(xf1(1:NFFT/2)));
% stem(k,2*abs(xf1(1:NFFT)))  
title('Espectro de Amplitud  de la señal')
xlabel('Frecuencia (Hz)')

waitfor(gcf)