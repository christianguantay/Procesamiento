% Filtro notch FLG
clear all;
close all;

%Colores
blue = [0, 0.4470, 0.7410];
orange = [0.8500, 0.3250, 0.0980];

%Leo la pista dada en el enunciado
[input,fs]=audioread('pista_01.wav');

%%Grafico de la interferencia
f1 = 210;
f2 = 375;
f3 = 720;

A1 = 0.05;
A2 = 0.03;
A3 = 0.02;

ts = 1/fs;

t = 0:ts:(length(input)/44100);
t = t(1:end-1);

noise = A1*cos(2*pi*f1*t) + A2*cos(2*pi*f2*t) + A3*cos(2*pi*f3*t) ;

figure,plot(t,noise,'Color', blue);
ylabel('Amplitud');
xlabel('t[s] (Tiempo)');
xlim([0 0.08]);
grid on;
str='Señal de interferencia';
title(str);
auxFFT  = abs(fft(noise,2^nextpow2(length(noise))));
auxFFT  = auxFFT/max(auxFFT); 
w       = 0:2*pi/(length(auxFFT)-1):2*pi;
figure,plot(w(1,1:floor(length(w)/2)),auxFFT(1,1:floor(length(w)/2)),'Color', blue);
ylabel('Amplitud');
xlabel('w(Frecuencia)');
ylim([0 1.2])
grid on;
str='Respuesta en frecuencia de la interferencia';
title(str);

%% Diseño los filtros notch

h1 = notch_flg_rect(0.0095);
h2 = notch_flg_rect(0.017);
h3 = notch_flg_rect(0.0326);

%% Creo  el filtro en cascada
h1f = dfilt.df2(h1,1);
h2f = dfilt.df2(h2,1);
h3f = dfilt.df2(h3,1);

hd = dfilt.cascade(h1f, h2f,h3f);

h12 = conv (h1,h2);
h123 = conv (h12,h3);

%% Filtrado de la señal con ruido
xconruido = input + noise';
y=conv(xconruido,h123); 
audiowrite('pista_01_noise.wav',xconruido,fs);
audiowrite('pista_01_cleaned.wav',y,fs);
y= y/max(y);

%% Graficos
%Macro usada para escalar los graficos a 1
MAG_MED = 15770;

X = fft(input,2^nextpow2(length(input)));
NOISE = fft(noise',2^nextpow2(length(input)));
XNOISE = fft(xconruido,2^nextpow2(length(input)));
Y = fft(y,2^nextpow2(length(y)));
w =(1:length(X)/2)*2*pi/length(X); % Vector de 0 a pi
H123 = fft(h123,2^nextpow2(length(input)));


figure(4);
plot(w/pi,abs(NOISE(1:length(w)))/MAG_MED, 'Color', orange);
grid on
hold on
plot(w/pi,abs(H123(1:length(w))),'Color', blue);
xlim([0 0.05])
title('Espectro de la interferencia')
xlabel('\omega [rad/s] (Frecuencia normalizada)')
ylabel('Amplitud')
str1 = 'Interferencia';
str2 = 'Respuesta en módulo del filtro '; 
legend(str1,str2,'Location','West')

%% Filtro en cascada. Graficos de amplitud y fase

figure(5)

plot(w/pi,abs(H123(1:length(w))),'Color', blue)
grid on
title('Filtro notch FLG diseñado con ventana rectangular')
xlabel('\omega [rad/s] (Frecuencia normalizada)')
ylabel('Amplitud')
xlim([0 0.1])
str = 'Respuesta en módulo del filtro total';
legend(str,'Location','NorthWest')

figure(6);
plot(w/pi,abs(XNOISE(1:length(w)))/MAG_MED, 'Color', orange);
grid on
hold on
plot(w/pi,abs(H123(1:length(w))),'Color', blue);
xlim([0 0.04])
title('Espectro de la señal original con ruido')
xlabel('\omega [rad/s] (Frecuencia normalizada)')
ylabel('Amplitud')
str1 = 'Señal original con ruido';
str2 = 'Respuesta en módulo del filtro '; 
legend(str1,str2,'Location','West')

figure(7);

plot(w/pi,abs(XNOISE(1:length(w)))/(max(abs(XNOISE(1:length(w))))),'Color', orange);
xlim([0 0.04])
ylim([0 1.2])
grid on
hold on

plot(w/pi,abs(Y(1:length(w)))/(max(abs(XNOISE(1:length(w))))),'Color', blue);
title('Espectro de la señal filtrada')
xlabel('\omega [rad/s] (Frecuencia normalizada)')
ylabel('Amplitud')
str1 = 'Señal original con ruido';
str2 = 'Señal filtrada '; 
legend(str1,str2,'Location','West')


function [hfilter] = notch_flg_rect(wc)
    %% Diseño del filtro notch
    % Creo este filtro como una combinacion de un filtro pasabajos y
    % pasaaltos
    
    %Desnormalizo wc
    wc = wc*pi;
    
    % Diseño del pasabajos
    wc_lpgen = wc - ((wc * 2.5)/100);
    wc_hpgen = wc + ((wc * 2.5)/100);
    
    
    ws_lpgen = wc_lpgen + ((wc_lpgen * 1)/100);      
    wp_lpgen = wc_lpgen - ((wc_lpgen * 1)/100); 
    
    %Redefino wc
    wc_lpgen = (ws_lpgen + wp_lpgen)/2;                

    delta_w_lpgen = ws_lpgen-wp_lpgen;                    

    % Diseño del pasaaltos
    ws_hpgen = wc_hpgen - ((wc_hpgen * 1)/100);      
    wp_hpgen = wc_hpgen + ((wc_hpgen * 1)/100);
    
    %Redefino wc
    wc_hpgen = (ws_hpgen + wp_hpgen)/2;  
    
    delta_w_hpgen = wp_hpgen-ws_hpgen;

    delta_wgen = min(delta_w_lpgen,delta_w_hpgen);

    N = ceil(4*pi/delta_wgen); % Orden del filtro. 
    if( mod(N,2)~=0 ) % Chequeo si es par
        N=N+1;
    end

    %Creo los filtros
    
    n = 0:N; 
    
    lpgen = sin(wc_lpgen*(n-N/2))./(pi*(n-N/2));   % sinc desplazada - Pasa bajos
    lpgen(N/2+1) = wc_lpgen/pi;
    
    hpgen = ((-1).^(n-N/2)).*sin((pi-wc_hpgen)*(n-N/2))./(pi*(n-N/2)); % Pasa altos
    hpgen((N/2+1)) = (pi-wc_hpgen)/pi;

    hfilter = hpgen+lpgen;
end