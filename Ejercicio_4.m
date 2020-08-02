clear all;
close all;
clc;

% Ejercicio 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cargo el filtro notch IIR y el ecualizador IIR
load notch_iir.mat;

load equalizador_iir.mat;

nfft = 4*2048;

% Se forma el numerador y denominador del sistema conjunto
num_neq = conv(num_eq,num_notch_iir);
den_neq = conv(den_eq,den_notch_iir);

[H_neq,w] = freqz(num_neq,den_neq,nfft);

W = [16,6,3];

num_neq_cuantizado = [];
den_neq_cuantizado = [];

for i=(1:length(W))
  num_neq_cuantizado(i,:) = round(num_neq*2^W(i))/(2^W(i));
  den_neq_cuantizado(i,:) = round(den_neq*2^W(i))/(2^W(i));
endfor

% Gráfico comparando los módulos

figure
hold on;
plot(w/pi,20*log10(abs(H_neq)));
for i=(1:length(W))
  hold on;
  [H,w] = freqz(num_neq_cuantizado(i,:),den_neq_cuantizado(i,:),nfft);
  plot(w/pi,20*log10(abs(H)));
endfor
grid on;
xlim([0 1]);
xlabel('\omega / \pi');
ylabel('Módulo [dB]');
legend('punto flotante',num2str(W(1)),num2str(W(2)),num2str(W(3)),'location','southeast');
title("Módulo del sistema ecualizador+notch IIR, para diferentes cuantizaciones");

% Gráfico de los diagramas de polos y zeros

% W = 16
figure
hold on;
polos_neq = roots(den_neq);
ceros_neq = roots(num_neq);
polos = roots(den_neq_cuantizado(1,:));
ceros = roots(num_neq_cuantizado(1,:));
zplane([ceros_neq, ceros],[polos_neq, polos]);
xlabel('Real');
ylabel('Imag');
title('Diagrama de polos y ceros del sistema cuantizado a 16 bits');

% W = 6
figure
hold on;
polos_neq = roots(den_neq);
ceros_neq = roots(num_neq);
polos = roots(den_neq_cuantizado(2,:));
ceros = roots(num_neq_cuantizado(2,:));
zplane([ceros_neq, ceros],[polos_neq, polos]);
xlabel('Real');
ylabel('Imag');
title('Diagrama de polos y ceros del sistema cuantizado a 6 bits');

% W = 4
figure
hold on;
polos_neq = roots(den_neq);
ceros_neq = roots(num_neq);
polos = roots(den_neq_cuantizado(3,:));
ceros = roots(num_neq_cuantizado(3,:));
zplane([ceros_neq, ceros],[polos_neq, polos]);
xlabel('Real');
ylabel('Imag');
title('Diagrama de polos y ceros del sistema cuantizado a 3 bits');

