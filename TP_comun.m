clear all;
close all;
clc;


% En este archivo se integran los diferentes filtros que se hicieron en los
% ejercicios, y se los prueba junto con las señales de audio.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Se importa el sistema electroacústico

load h_sys.mat;
fs = 44.1e3;
nfft = 4*2048;

k = 2;



% Ejercicio 2
% Se importa el filtro Notch FIR

load notch_fir.mat;

% Se importan los filtros Notch IIR

load notch_iir.mat;

% Ejercicio 3
% Se importa el equalizador FIR

load equalizador_fir.mat;

% Se importa el equalizador IIR
% Se importa en forma de coeficientes del numerador y denominador
load equalizador_iir.mat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Se lee el archivo de pista
[x,fs]=audioread('CANCIONES/pista_01.wav');
x = x';

% Se leen las interferencias
%load interferencia.mat;

% Se le aplica el sistema electroacústico

x_1 = filter(h,1,x);

% Se suman las interferencias a la señal
f1 = 210;
f2 = 375;
f3 = 720;

A1 = 0.05;
A2 = 0.03;
A3 = 0.02;

ts = 1/fs;

t = 0:ts:(length(x)/fs);
t = t(1:end-1);

noise = A1*cos(2*pi*f1*t) + A2*cos(2*pi*f2*t) + A3*cos(2*pi*f3*t);

x_2 = x_1 + noise;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(k==1)
% Primer caso: filtro notch FIR + ecualizador FIR

% Se aplica el filtro notch FIR

x_31 = filter(h123,1,x_2);

% Se aplica el ecualizador FIR

x_41 = filter(h_eq,1,x_31);

% El sistema total queda de la siguiente forma:

h_tot_1 = conv(h,h123);
h_tot_1 = conv(h_tot_1,h_eq);

[H_tot,W_tot] = freqz(h_tot_1,1,nfft);

figure
hold on;
plot(W_tot/pi,20*log10(abs(H_tot)),'b-','linewidth',1);
hold on;
line([0 1], [2 2], "linestyle","-","color","k","linewidth",1);
hold on;
line([0 1], [-2 -2], "linestyle","-","color","k","linewidth",1);
grid on;
legend('Sistema total','Límite +2 dB','Límite -2 dB')
xlabel('\omega / \pi (Frecuencia normalizada)');
ylabel('Amplitud [dB]');
title('Sistema total');
%
%% Se reproduce la señal en todas sus etapas
%sound(x,fs);    % Señal original
%sound(x_1,fs);  % Señal luego del sistema electroacústico
%sound(x_2,fs);  % Señal luego del sistema electroacústico más ruido
%sound(x_31,fs); % Señal luego del notch
%sound(x_41,fs); % Señal luego del filtro ecualizador
endif
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Segundo caso: filtro notch IIR + ecualizador IIR
if(k==2)
% Se aplica el primer filtro Notch IIR

x_321 = filter(numerador_f1,denominador_f1,x_2);
x_322 = filter(numerador_f2,denominador_f2,x_321);
x_32 = filter(numerador_f3,denominador_f3,x_322);

%Se reproduce la señal en todas sus etapas
%sound(x,fs);    % Señal original
%sound(x_1,fs);  % Señal luego del sistema electroacústico
%sound(x_2,fs);  % Señal luego del sistema electroacústico más ruido
%sound(x_32,fs); % Señal luego del notch

y = fft(noise,nfft/4);


w = linspace(0,2*pi,nfft/4);

figure
hold on;
plot(w/pi,abs(y));
xlim([0 1]);
grid on;

endif