clear all;
close all;
clc;

% Ejercicio 3.2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load("h_sys.mat");
fs = 44.1e3;
nfft = 4*2048;

% Diagrama de polos y ceros del sistema electroacústico
figure
hold on;
zplane(h,1);
xlabel("Real");
ylabel("Imag");
title("Diagrama de polos y zeros del sistema electroacústico");

% Cálculo de los ceros del sistema H

ceros_H = roots(fliplr(h));

M = length(ceros_H);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diagrama de polos y ceros del pasa todo y del fase mínima
% Se buscan los polos y ceros de cada uno, en función de los del sistema
% original

ceros_H_min = [];
polos_H_min = zeros(1,M);

ceros_H_ap = [];
polos_H_ap = [];


% Se separan los ceros de cada uno de los filtros: pasa todo y de fase mínima
for i=(1:M)
  Z0 = ceros_H(i);
  if abs(Z0)>1
    ceros_H_min = [ceros_H_min inv(conj(Z0))];
    ceros_H_ap = [ceros_H_ap Z0];
    polos_H_ap = [polos_H_ap inv(conj(Z0))];
  else
    ceros_H_min = [Z0 ceros_H_min];
  endif
endfor

% Filtro pasa todo

num_ap = prod((-1)*(1./ceros_H_ap))*poly(ceros_H_ap);
den_ap = poly(polos_H_ap);

figure
hold on;
zplane(num_ap,den_ap);
xlabel("Real");
ylabel("Imag");
title("Diagrama de polos y zeros del sistema pasa todo");

% Filtro fase minima

num_min = poly(ceros_H_min)/( prod((-1)*(1./ceros_H_ap))/h(114));
den_min = poly(polos_H_min);

figure
hold on;
zplane(num_min,den_min);
title("Diagrama de polos y zeros del sistema fase mínima");

% Sistema ecualizador que es el inverso del de fase mínima

num_eq = den_min;
den_eq = num_min;

figure
hold on;
zplane(num_eq,den_eq);
title("Diagrama de polos y zeros del sistema ecualizador");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Respuesta en frecuencia de todos los sistemas
[H,W] = freqz(h,1,nfft);

H_eq = freqz(num_eq,den_eq,W');

% Sistema total
num_tot = conv(h,num_eq);
den_tot = conv(1,den_eq);


H_tot = freqz(num_tot,den_tot,W');

% Gráfico de comparando el módulo de los sistemas
figure
hold on;
plot(W/pi,20*log10(abs(H)));
hold on;
plot(W/pi,20*log10(abs(H_eq)));
hold on;
plot(W/pi,20*log10(abs(H_tot)));
grid on;
xlim([0 1]);
xlabel('\omega / \pi');
ylabel('Módulo [dB]');
title("Módulo de los diferentes sistemas");
legend("Sistema electroacústico","Ecualizador","Sistema total");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Retardo de grupo de los sistemas

[G, W] = grpdelay(h,1,nfft); %Retardo de grupo del sistema electroacúsitco

[G_eq, W_eq] = grpdelay(num_eq,den_eq,nfft); %Retardo de grupo del sistema ecualizador

[G_tot, W_tot] = grpdelay(num_tot,den_tot,nfft); %Retardo de grupo del sistema total

% Gráfico comparando los retardos de grupo

figure
hold on;
plot(W/pi,G);
hold on;
plot(W_eq/pi,G_eq);
hold on;
plot(W_eq/pi,G_tot);
grid on;
xlim([0 1]);
xlabel('\omega / \pi');
ylabel('Muestras');
title("Retardo de grupo de los diferentes sistemas");
legend("Sistema electroacústico","Ecualizador","Sistema total");

% Se exportan los datos del numerador y denominador del filtro

save equalizador_iir.mat num_eq den_eq;
