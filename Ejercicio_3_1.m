clear all;
close all;
clc;


% Ejercicio 3.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load h_sys.mat;
fs = 44.1e3;

nfft = 4*2048;

[H,W] = freqz(h,1,nfft);

% Se define el sistema inverso

[H_inv,W_inv] = freqz(1,h,nfft);

delta_H_tot = 0.2;

figure
hold on;
plot(W_inv/pi,abs(H_inv),'linewidth',1,'b-');
%hold on;
%plot(W_inv/pi,abs(H_inv)*(1+delta_H_tot),'linewidth',1,'r-');
%hold on;
%plot(W_inv/pi,abs(H_inv)*(1-delta_H_tot),'linewidth',1,'r-');
xlim([0 1]);
xlabel('\omega / \pi (Frecuencia normalizada)');
ylabel('Amplitud');
title("Respuesta en frecuencia del sistema inverso");
grid on;

% A partir del gráfico veo los límites para cada banda

% Banda 1
A1_sup = 1.37838;
%A1_sup = 1.5;
A1_inf = 0.925499;
%A1_inf = 1.1;

n_1_inicio = 1;  % Inicio de la primer banda
%n_1_fin = ceil(0.0232 * nfft);  % Fin de la primera banda
n_1_fin = ceil(0.035 * nfft);

% Banda 2
A2_sup = 0.227298;
A2_inf = 0.152176;

%n_2_inicio = ceil(0.0862670 * nfft);
n_2_inicio = ceil(0.048 * nfft);
%n_2_fin = ceil(0.0963561 * nfft);
n_2_fin = ceil(0.122 * nfft);

% Banda 3
A3_sup = 4.45815;
A3_inf = 3.15189;

%n_3_inicio = ceil(0.156569 * nfft);
n_3_inicio = ceil(0.178 * nfft);
%n_3_fin = ceil(0.425339 * nfft);
n_3_fin = ceil(0.385*nfft);

% Banda 4
A4_sup = 0.68;
A4_inf = 0.474943;

%n_4_inicio = ceil(0.49 * nfft);
n_4_inicio = ceil(0.492 * nfft);
%n_4_fin = ceil(0.69 * nfft);
n_4_fin = ceil(0.72 * nfft);

% Banda 5
A5_sup = 1.95604;
A5_inf = 1.37825;

%n_5_inicio = ceil(0.76 * nfft);
n_5_inicio = ceil(0.75 * nfft);
%n_5_fin = ceil(0.782288 * nfft);
n_5_fin = ceil(0.79 * nfft);

% Banda 6
A6_sup = 5;
A6_inf = 3.52376;

%n_6_inicio = ceil(0.965 * nfft);
n_6_inicio = ceil(0.95*nfft);
n_6_fin = ceil(1 * nfft);

figure
hold on;
plot(W_inv/pi,abs(H_inv),'linewidth',1,'b-');
hold on;
plot(W_inv/pi,abs(H_inv)*(1+delta_H_tot),'linewidth',1,'r-');
hold on;
plot(W_inv/pi,abs(H_inv)*(1-delta_H_tot),'linewidth',1,'r-');


n_inicio = [n_1_inicio, n_2_inicio, n_3_inicio, n_4_inicio, n_5_inicio, n_6_inicio];
n_fin = [n_1_fin, n_2_fin, n_3_fin, n_4_fin, n_5_fin, n_6_fin];
A_sup = [A1_sup, A2_sup, A3_sup, A4_sup, A5_sup, A6_sup];
A_inf = [A1_inf, A2_inf, A3_inf, A4_inf, A5_inf, A6_inf];

for i=(1:6)
  hold on;
  line([W_inv(n_inicio(i))/pi W_inv(n_fin(i))/pi], [A_sup(i) A_sup(i)], "linestyle","-","color","k","linewidth",1);
  hold on;
  line([W_inv(n_inicio(i))/pi W_inv(n_fin(i))/pi], [A_inf(i) A_inf(i)], "linestyle","-","color","k","linewidth",1);
  hold on;
  x =  W_inv(n_inicio(i))/pi;
  plot([x,x],[A_sup(i), 0],'k-','linewidth',1);
  hold on;
  x =  W_inv(n_fin(i))/pi;
  plot([x,x],[A_sup(i), 0],'k-','linewidth',1);
endfor

xlim([0 1]);
xlabel('\omega / \pi (Frecuencia normalizada)');
ylabel('Amplitud');
title("Bandas identificadas en el sistema inverso");
grid on;
legend('sistema inverso','tolerancia máxima permitida');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diseño del multibanda
delta_bandas = (A_sup - A_inf)/2;

delta = min(delta_bandas);  % Diseño para una ventana con este delta

% Las banda son contiguas luego elijo la ventana considerando esto
A_med = (A_sup + A_inf)/2;
A_med(2) = 0.22;
A_med(3) = 4.3;
A_med(5) = 1.9;
A_med(6) = 3.6;
max_dist_bandas = 0;

min_delta_w = 1000;

for i=(1:5)
  k = A_med(i)-A_med(i+1);
  w = W_inv(n_fin(i+1)) - W_inv(n_inicio(i));
  if k>max_dist_bandas
    max_dist_bandas = k;
  endif
  if w < min_delta_w
    min_delta_w = w;
  endif
endfor

delta_ventana = delta / max_dist_bandas;  % A partir de este valor elijo la ventana

% Utilizo la ventana de Hann

% Busco la distancia entre bandas mínima para elegir el orden de la ventana


M = ceil(8 * pi / min_delta_w);

M = 112;

if mod(M,2)!=0
  M = M+1;
endif
vent = window(@hamming, M+1)';

% Defino el multibanda ideal

n = [0:1:M];

h_mb = zeros(1,M+1);

A = [A_med 0];

n_corte = zeros(1,6);
for i=(1:5)
  n_corte(i) = (n_fin(i) + n_inicio(i+1))/2;
endfor
n_corte(6) = nfft;

w_corte = W_inv(floor(n_corte));

for i=(1:6)
  h_mb = h_mb + (A(i) - A(i+1)) .* sin(w_corte(i) .* (n.- M/2))./(pi * (n .- M/2));
endfor

h_mb(M/2 + 1) = 0;
for i=(1:6)
  h_mb(M/2 + 1) = h_mb(M/2 + 1) + (A(i) - A(i+1))*w_corte(i)/pi;
endfor

h_eq = h_mb.*vent;

[H_eq, W_eq] = freqz(h_eq,1,nfft);

figure
hold on;
plot(W_eq/pi, abs(H_eq));
%hold on;
%plot(W_inv/pi,abs(H_inv)*(1+delta_H_tot),'linewidth',1,'r-');
%hold on;
%plot(W_inv/pi,abs(H_inv)*(1-delta_H_tot),'linewidth',1,'r-');
hold on;
%plot(W_inv/pi,abs(H_inv),'linewidth',1,'r-');
legend('multibanda diseñado','Sistema inverso');
for i=(1:6)
  hold on;
  line([W_inv(n_inicio(i))/pi W_inv(n_fin(i))/pi], [A_med(i)+delta A_med(i)+delta], "linestyle","-","color","k","linewidth",1);
  hold on;
  line([W_inv(n_inicio(i))/pi W_inv(n_fin(i))/pi], [A_med(i)-delta A_med(i)-delta], "linestyle","-","color","k","linewidth",1);
  hold on;
  x =  W_inv(n_inicio(i))/pi;
  plot([x,x],[A_med(i)+delta, 0],'k-','linewidth',1);
  hold on;
  x =  W_inv(n_fin(i))/pi;
  plot([x,x],[A_med(i)+delta, 0],'k-','linewidth',1);
endfor
xlim([0,1]);
grid on;
xlabel('\omega / \pi (Frecuencia normalizada)');
ylabel('Amplitud');
title(['Filtro multibanda diseñado, de orden ',num2str(M)]);

% Se fueron ajustando los valores tanto de las bandas como de las frecuencias 
% de corte de cada una de las bandas, de forma de obtener resultados
% aceptables.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Resultado final

%H_tot = H.*H_eq;

[H_tot, W_tot] = freqz(conv(h,h_eq),1,nfft);

%Comparación de los resultados

figure
hold on;
plot(W_inv/pi,20*log10(abs(H_tot)),'b-','linewidth',1);
hold on;
line([0 1], [2 2], "linestyle","-","color","k","linewidth",1);
hold on;
line([0 1], [-2 -2], "linestyle","-","color","k","linewidth",1);
hold on;
plot(W_eq/pi,20*log10(abs(H_eq)),'r-','linewidth',1);
hold on;
plot(W/pi,20*log10(abs(H)),'g-','linewidth',1);
grid on;
xlim([0 1]);
grid on;
xlabel('\omega / \pi (Frecuencia normalizada)');
ylabel('Módulo [dB]');
title("Comparación del módulo de los sistemas");
set(gca, 'ytick', -15:1:15);
legend("Sistema total","Límite +2 dB","Límite -2 dB","Sistema equalizador", "Sistema Electroacústico");

% Comparación de los retardos de fase

% Sistema electroacústico
fase_h = atan2(imag(H),real(H));

for i=(2:nfft)
  if abs(fase_h(i)-fase_h(i-1))>=1.9*pi
    fase_h = [fase_h(1:i-1) ; (fase_h(i:end) .- 2*pi)];
  endif
endfor

retardo_fase_h = fase_h./W;

% Sistema ecualizador
retardo_fase_h_eq = ones(nfft,1).*(M/2);  % Se calcula así porque es FLG tipo 1

% Sistema total
retardo_fase_h_tot = retardo_fase_h + retardo_fase_h_eq;  % Se suman porque
% están en cascada

figure
hold on;
plot(W/pi,retardo_fase_h_tot,'b-','linewidth',1);
hold on;
plot(W/pi,retardo_fase_h_eq,'r-','linewidth',1);
hold on;
plot(W/pi,retardo_fase_h,'g-','linewidth',1);
grid on;
xlabel('\omega / \pi (Frecuencia normalizada)');
ylabel('Muestras');
title("Comparación de los retardos de fase");
legend("Sistema total","Sistema equalizador", "Sistema Electroacústico");

% Se exporta el sistema equalizador FIR 

save equalizador_fir.mat h_eq;
