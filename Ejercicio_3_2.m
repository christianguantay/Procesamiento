clear all;
close all;
clc;

% Ejercicio 3.2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load("h_sys.mat");
fs = 44.1e3;

% Diagrama de polos y ceros del sistema electroacústico
figure
hold on;
zplane(h,1);
xlabel("Real");
ylabel("Imag");
title("Diagrama de polos y zeros del sistema electroacústico");

% Cálculo de los ceros del sistema H

ceros_H = roots(h);

M = length(ceros_H);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diagrama de polos y ceros del pasa todo y del fase mínima
% Se buscan los polos y ceros de cada uno, en función de los del sistema
% original

ceros_H_min = [];
polos_H_min = zeros(1,M);

ceros_H_ap = [];
polos_H_ap = [];

for i=(1:M)
  Z0 = ceros_H(i);
  if abs(Z0)>1
    ceros_H_min = [ceros_H_min inv(conj(Z0))];
    ceros_H_ap = [ceros_H_ap Z0];
    polos_H_ap = [polos_H_ap inv(conj(Z0))];
  else
    ceros_H_min = [ceros_H_min Z0];
  endif
endfor

% Filtro pasa todo

num_ap = poly(ceros_H_ap);
den_ap = poly(polos_H_ap);

figure
hold on;
zplane(num_ap,den_ap);
xlabel("Real");
ylabel("Imag");
title("Diagrama de polos y zeros del sistema pasa todo");

% Filtro fase minima

num_min = poly(ceros_H_min);
den_min = poly(polos_H_min);

figure
hold on;
plot(real(ceros_H_min),imag(ceros_H_min),'bo');
plot(real(polos_H_min),imag(polos_H_min),'bx');
hold on;
t = linspace(0,2*pi,100)'; 
circsx = 1.*cos(t); 
circsy = 1.*sin(t); 
plot(circsx,circsy); 
grid on;
xlabel("Real");
ylabel("Imag");
title("Diagrama de polos y zeros del sistema fase mínima");

% Sistema ecualizador que es el inverso del de fase mínima

ceros_H_eq = polos_H_min;
polos_H_eq = ceros_H_min;

figure
hold on;
plot(real(ceros_H_eq),imag(ceros_H_eq),'bo');
plot(real(polos_H_eq),imag(polos_H_eq),'bx');
hold on;
t = linspace(0,2*pi,100)'; 
circsx = 1.*cos(t); 
circsy = 1.*sin(t); 
plot(circsx,circsy); 
grid on;
xlabel("Real");
ylabel("Imag");
title("Diagrama de polos y zeros del sistema ecualizador");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Respuesta en frecuencia de todos los sistemas

[H,W] = freqz(h,1,512,'whole');

num_eq = conv(num_ap,1) * prod(polos_H_ap);
den_eq = conv(den_ap,h);
% Como el cálculo de los coeficientes del de fase mínima daba problemas
% usando la función poly, se expresa el numerador y el numerador del
% sistema ecualizador en función del sistema original y el pasa todo

[H_eq,W_eq] = freqz(num_eq ,den_eq);

% Sistema total
num_tot = conv(h,num_eq);
den_tot = conv(1,den_eq);

[H_tot, W_tot] = freqz(num_tot,den_tot);

figure
hold on;
plot(W/pi,20*log10(abs(H)));
hold on;
plot(W_eq/pi,20*log10(abs(H_eq)));
hold on;
plot(W_tot/pi,20*log10(abs(H_tot)));
grid on;
xlim([0 1]);
xlabel('\omega / \pi');
ylabel('Módulo [dB]');
title("Módulo de los diferentes sistemas");
legend("Sistema electroacústico","Ecualizador","Sistema total");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Retardo de grupo de los sistemas

[G, W] = grpdelay(h); %Retardo de grupo del sistema electroacúsitco

[G_eq, W_eq] = grpdelay(num_eq,den_eq); %Retardo de grupo del sistema ecualizador

[G_tot, W_tot] = grpdelay(num_tot,den_tot); %Retardo de grupo del sistema total

%G = G./fs;
%G_eq = G_eq ./fs;
%G_tot = G_tot./fs;

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
ylabel('Retardo de grupo [ms]');
title("Retardo de grupo de los diferentes sistemas");
legend("Sistema electroacústico","Ecualizador","Sistema total");

save equalizador_iir.mat num_eq den_eq;



