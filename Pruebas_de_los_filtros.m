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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Esta constante sirve solamente para compilar distintas partes de este archivo
% k = 1 ==> Ecualizador FIR
% k = 2 ==> Ecualizador IIR
% k = 3 ==> Comparación de los ecualizadores
k = 1; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ejercicio 2
% Se importa el filtro Notch FIR

load notch_fir.mat;

% Se importan los filtros Notch IIR

load notch_iir.mat;

% Ejercicio 3
% Se importa el equalizador FIR

load equalizador_fir.mat;

% Se importa el equalizador IIR
load equalizador_iir.mat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Se lee el archivo de pista
[x,fs]=audioread('CANCIONES/pista_06.wav');

x = x';

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
L = length(x);
t = (0:L-1) * ts;

noise = A1*cos(2*pi*f1*t) + A2*cos(2*pi*f2*t) + A3*cos(2*pi*f3*t);

x_2 = x_1 + noise;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Primer caso: electroacústico + ecualizador FIR
if(k==1)

% Se aplica el ecualizador FIR

x_41 = filter(h_eq,1,x_1);

% Se reproduce la señal en todas sus etapas
sound(x,fs);    % Señal original
sound(x_1,fs);  % Señal luego del sistema electroacústico
sound(x_41,fs); % Señal luego del filtro ecualizador


endif
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Segundo caso: ecualizador IIR
if(k==2)

% Se aplica el ecualizador IIR

x_42 = filter(num_eq,den_eq,x_1);

% Se reproduce la señal en todas sus etapas
sound(x,fs);    % Señal original
sound(x_1,fs);  % Señal luego del sistema electroacústico
sound(x_42,fs); % Señal luego del filtro ecualizador

endif
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tercer caso: Se compara a los ecualizadores
if(k==3)
% Se aplica el ecualizador IIR y el FIR y se reproducen ambas señales

x_431 = filter(h_eq,1,x_1);
x_432 = filter(num_eq,den_eq,x_1);

%Se reproduce la señal en todas sus etapas
sound(x,fs);    % Señal original
sound(x_431,fs); % Señal luego del filtro ecualizador FIR
sound(x_432, fs); % Señal luego del filtro ecualizador IIR

endif