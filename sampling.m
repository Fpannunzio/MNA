clear all ; close all; clc

Fs = 25;                      % Sampling rate
T = 1/Fs;                     % Sample time
L = 4000;                     % Length of signal
t = (0:L-1)*T;                % Time vector

A1=1;
A2= 3;


NP=2;

y = A1*sin(40*pi*t) + A2*sin(10*pi*t).^2;

y = A1*sin(30*pi*t) + A2*sin(20*pi*t);

%%%%%%%%%%%%%%%%%%%%%%
NFFT = 2^nextpow2(NP*L); % Next power of 2 from length of y
Y = fft(y,NFFT);
f = Fs/2*linspace(0,1,NFFT/2);
% Plot single-sided amplitude spectrum.
stem(f,2*abs(Y(1:NFFT/2))/L)   
title('Espectro de Amplitud  de la señal')
xlabel('Frecuencia (Hz)')

waitfor(gcf)