clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[numerador,denominador] = filtro_notch_iir_butterworth(f1,r1,r2,delta_p,delta_s,fs)
%function[z,p,g] = filtro_notch_iir_butterworth(f1,r1,r2,delta_p,delta_s,fs)
% Recibe "f", la frecuencia que se quiere eliminar, r es un coeficiente 
% que define a qué distancia se encuentran las frecuencias de transición
% a cada lado de la frecuencia "f", "delta_p" es la tolerancia en la banda de
% paso y "delta_s" es la tolerancia en la banda eliminada. El parámetro "fs" es
% la frecuencia de muestreo para la aplicación de la transformación bilineal.
% Devuelve los coeficientes del numerador y el denominador del filtro Notch IIR 
% discreto a partir de un filtro pasa bajos continuo Butterworth.
    
  omega_notch_p_1 = (f1-f1*r1) * 2 * pi;
  omega_notch_s_1 = (f1-f1*r2) * 2 * pi;

  omega_notch_s_2 = (f1+f1*r2) * 2 * pi;
  omega_notch_p_2 = (f1+f1*r1) * 2 * pi;
  
  % Con esto defino las especificaciones del pasa bajos

  omega_l = omega_notch_p_1;
  omega_h = omega_notch_p_2;
  omega_p = 1;

  % Calculo los omega_s_1 y omega_s_2 y defino omega_s como el mínimo de ellos
 
  omega_s_1 = omega_notch_s_1 * (omega_h - omega_l)/(omega_l * omega_h - omega_notch_s_1^2);
  omega_s_2 = omega_notch_s_2 * (omega_h - omega_l)/(omega_l * omega_h - omega_notch_s_2^2);

  omega_s = min(abs(omega_s_1),abs(omega_s_2));

  % Calculo los coeficientes "d" y "k" para el pasa bajos

  d = sqrt(((1-delta_p)^(-2) -1)/(delta_s^(-2) - 1));
  k = omega_p/omega_s;

  N = ceil(log(1/d)/log(1/k));  % Orden del pasa bajos
  omega_0 = (omega_p * ((1-delta_p)^(-2) -1)^(-1/(2*N)) + omega_s * (delta_s^(-2) - 1)^(-1/(2*N)))/2;

  % Butterworth pasa bajos IIR continuo
  
  [z,p,g] = butter(N,omega_0,"low","s");
% Con esto defino los polos y ceros del notch IIR

  g_notch = g / prod((-p));
  
  cero_del_notch = [j * sqrt(omega_l * omega_h), -j * sqrt(omega_l * omega_h)];

  z_notch = [];
  p_notch = [];

  for i=(1:N)
    z_notch = [z_notch cero_del_notch];
    %polo_del_notch = [(p(i)^(-1) + sqrt(p(i)^(-2)*(omega_h - omega_l)^(2) -4 * omega_h * omega_l))/2, (p(i)^(-1) - sqrt(p(i)^(-2)*(omega_h - omega_l)^(2) -4 * omega_h * omega_l))/2];
    polo_del_notch = [((-omega_h + omega_l) + sqrt((omega_h - omega_l)^2 -4 * p(i)^2 * omega_h * omega_l))/(-2 * p(i)) , ((-omega_h + omega_l) - sqrt((omega_h - omega_l)^2 -4 * p(i)^2 * omega_h * omega_l))/(-2 * p(i))];
    p_notch = [p_notch polo_del_notch];
  endfor

  % Transformo usando la transformada bilineal

  [zb_notch, za_notch] = bilinear(z_notch,p_notch,g_notch,1/fs);
  numerador = zb_notch;
  denominador = za_notch;
  %[z,p,g] = bilinear(z_notch,p_notch,g_notch,1/fs);
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filtro Notch IIR
fs = 44.1e3;
% Frecuencias que se quieren eliminar
f1 = 210;
f2 = 375;
f3 = 720;
% Puntos para el cálculo de las respuestas en frecuencia
nfft = 8*2048;
% Parámetros de diseño de los filtros
r1= 0.1;
r2 = 0.0025;
delta_p = 0.02;
delta_s = 0.1;



% Se calculan los filtros IIR Notch discretos

%[z_1,p_1,g_1] = filtro_notch_iir_butterworth(f1,r1,r2,delta_p,delta_s,fs);
%[z_2,p_2,g_2] = filtro_notch_iir_butterworth(f2,r1,r2,delta_p,delta_s,fs);
%[z_3,p_3,g_3] = filtro_notch_iir_butterworth(f2,r1,r2,delta_p,delta_s,fs);

[numerador_f1, denominador_f1] = filtro_notch_iir_butterworth(f1,r1,r2,delta_p,delta_s,fs);
[numerador_f2, denominador_f2] = filtro_notch_iir_butterworth(f2,r1,r2,delta_p,delta_s,fs);
[numerador_f3, denominador_f3] = filtro_notch_iir_butterworth(f3,r1,r2,delta_p,delta_s,fs);

% Gráfico de los 3 filtros juntos

[H1,W1] = freqz(numerador_f1,denominador_f1,nfft);
[H2,W2] = freqz(numerador_f2,denominador_f2,nfft);
[H3,W3] = freqz(numerador_f3,denominador_f3,nfft);
N = length(numerador_f1);
figure
hold on;
plot(W1/pi,abs(H1),'b-','linewidth',1);
hold on;
plot(W2/pi,abs(H2),'r-','linewidth',1);
hold on;
plot(W3/pi,abs(H3),'g-','linewidth',1);
grid on;
xlim([0 0.06]);
ylim([0 1.1]);
xlabel('\omega / \pi');
ylabel('Amplitud');
leyenda1 = sprintf('f1 = %i',f1);
leyenda2 = sprintf('f2 = %i',f2);
leyenda3 = sprintf('f3 = %i',f3);
legend(leyenda1,leyenda2,leyenda3,'location','southeast');
titulo = sprintf('Respuesta en frecuencia de los filtros notch IIR de orden %i',N);
title(titulo);


% Se exportan los filtros IIR para luego usarlos

save notch_iir.mat numerador_f1 denominador_f1 numerador_f2 denominador_f2 numerador_f3 denominador_f3;

