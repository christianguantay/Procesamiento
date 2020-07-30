% Filtro notch FLG
clear all;
close all;

%Leo la pista dada en el enunciado
[x,fs]=audioread('pista_01.wav');

%%Grafico de la interferencia
f1 = 210;
f2 = 375;
f3 = 720;

A1 = 0.05;
A2 = 0.03;
A3 = 0.02;

ts = 1/fs;

t = 0:ts:(length(x)/44100);
t = t(1:end-1);

nfft = 4*2048;

noise = A1*cos(2*pi*f1*t) + A2*cos(2*pi*f2*t) + A3*cos(2*pi*f3*t) ;
[h,w]=freqz(noise,nfft);

figure, plot(w,abs(h),'b');

figure,plot(t,noise,'r');
xlim([0 0.08]);
grid on;
str='0.08s de interferencia';
title(str);
auxFFT  = abs(fft(noise,2^nextpow2(length(noise))));
w       = 0:2*pi/(length(auxFFT)-1):2*pi;
figure,plot(w(1,1:floor(length(w)/2)),auxFFT(1,1:floor(length(w)/2)),'r');
grid on;
str='FFT de la interferencia';
title(str);

%% Dise�o del filtro notch 1
% Creo este filtro como una combinacion de un filtro pasabajos y pasaalto

% Dise�o del pasabajos 1

ws_lp = 0.0094*pi;      
wp_lp = 0.009*pi; %0.0092
wc_lp = (ws_lp+wp_lp)/2;                

delta_w_lp=ws_lp-wp_lp;                    

% Dise�o del pasaaltos 1
ws_hp = 0.0096*pi;      
wp_hp = 0.01*pi; %0.0098
wc_hp = (ws_hp+wp_hp)/2;  

delta_w_hp=wp_hp-ws_hp;

delta_w=min(delta_w_lp,delta_w_hp);

N= ceil(4*pi/delta_w); % Orden del filtro.
if( mod(N,2)~=0 ) % Chequeo si es par
    N=N+1;
end
%Creo los filtros

n = 0:N; 

sinc_desp=sin(wc_lp*(n-N/2))./(pi*(n-N/2));   % sinc desplazada - Pasa bajos
sinc_desp(N/2+1)=wc_lp/pi;
lp = sinc_desp;                               % para mayor claridad en el codigo

q = pi-wc_hp;
hp=((-1).^(n-N/2)).*sin(q*(n-N/2))./(pi*(n-N/2)); % Pasa altos
hp((N/2+1))=q/pi;

h=lp+hp; 

%% Dise�o del filtro notch 2
%Creo este filtro como una combinacion de un filtro pasabajos y pasaalto

%Dise�o del pasabajos 2

ws_lp = 0.0168*pi;      
wp_lp = 0.0164*pi; 
wc_lp = (ws_lp+wp_lp)/2;                

delta_w_lp=ws_lp-wp_lp;                    

% Dise�o del pasaaltos 2
ws_hp = 0.0172*pi;      
wp_hp = 0.0176*pi; 
wc_hp = (ws_hp+wp_hp)/2;  

delta_w_hp=wp_hp-ws_hp;

delta_w=min(delta_w_lp,delta_w_hp);

N= ceil(4*pi/delta_w); % Orden del filtro.
if( mod(N,2)~=0 ) % Chequeo si es par
    N=N+1;
end

%Creo los filtros

n = 0:N; 

sinc_desp=sin(wc_lp*(n-N/2))./(pi*(n-N/2));   % sinc desplazada - Pasa bajos
sinc_desp(N/2+1)=wc_lp/pi;
lp2 = sinc_desp;                               % para mayor claridad en el codigo

q = pi-wc_hp;
hp2=((-1).^(n-N/2)).*sin(q*(n-N/2))./(pi*(n-N/2)); % Pasa altos
hp2((N/2+1))=q/pi;

h2=lp2+hp2;

%% Dise�o del filtro notch
% Creo este filtro como una combinacion de un filtro pasabajos y pasaalto

% Dise�o del pasabajos 3

ws_lp = 0.0325*pi;      
wp_lp = 0.0305*pi; 
wc_lp = (ws_lp+wp_lp)/2;                

delta_w_lp=ws_lp-wp_lp;                    

% Dise�o del pasaaltos 3
ws_hp = 0.0335*pi;      
wp_hp = 0.0355*pi; 
wc_hp = (ws_hp+wp_hp)/2;  
delta_w_hp=wp_hp-ws_hp;

delta_w=min(delta_w_lp,delta_w_hp);

N= ceil(4*pi/delta_w); % Orden del filtro. 
if( mod(N,2)~=0 ) % Chequeo si es par
    N=N+1;
end

%Creo los filtros
n = 0:N; 

sinc_desp=sin(wc_lp*(n-N/2))./(pi*(n-N/2));   % sinc desplazada - Pasa bajos
sinc_desp(N/2+1)=wc_lp/pi;
lp3 = sinc_desp;                               % para mayor claridad en el codigo

q = pi-wc_hp;
hp3=((-1).^(n-N/2)).*sin(q*(n-N/2))./(pi*(n-N/2)); % Pasa altos
hp3((N/2+1))=q/pi;

h3=hp3+lp3; 

%Creo  el filtro en cascada
%h1f = dfilt.df2(h,1);
%h2f = dfilt.df2(h2,1);
%h3f = dfilt.df2(h3,1);

%hd = dfilt.cascade(h1f, h2f,h3f);

h12 = conv (h,h2);
h123 = conv (h12,h3);
%islinphase(h123)

%% Filtrado de la se�al con ruido
xconruido = x + noise';

y=conv(xconruido,h123); 
audiowrite('pista_01_noise.wav',xconruido,fs);
audiowrite('pista_01_cleaned.wav',y,fs);
y= y/max(y);

%% Graficos

X = fft(x,2^nextpow2(length(x)));

XNOISE = fft(xconruido,2^nextpow2(length(x)));

Y = fft(y,2^nextpow2(length(y)));

w =(1:length(X)/2)*2*pi/length(X); % Vector de 0 a pi

%% Filtro notch 1

LP = fft(lp,2^nextpow2(length(x)));

HP = fft(hp,2^nextpow2(length(x)));

H = fft(h,2^nextpow2(length(x)));

%Filtro notch 2

LP2 = fft(lp2,2^nextpow2(length(x)));

HP2 = fft(hp2,2^nextpow2(length(x)));

H2 = fft(h2,2^nextpow2(length(x)));

%Filtro notch 3

LP3 = fft(lp3,2^nextpow2(length(x)));

HP3 = fft(hp3,2^nextpow2(length(x)));

H3 = fft(h3,2^nextpow2(length(x)));

%Filtro en cascada
H123 = fft(h123,2^nextpow2(length(x)));


%% Filtro en cascada. Graficos de amplitud y fase
figure,plot(w/pi,abs(H123(1:length(w))));
grid on;
title('Filtro pasabajos dise�ado con ventana rectangular');
xlabel('\omega (rad/s)');
ylabel('Amplitud');
xlim([0 0.1]);
str = 'Respuesta en m�dulo del filtro total';
legend(str,'Location','NorthWest');


%figure, plot(w/pi, angle(H123(1:length(w)))), title('Phase plot')

%Se�al filtrada y se�al con ruido

figure()

subplot(2,1,1);
plot(w/pi,abs(XNOISE(1:length(w))));
grid on;
hold on;
plot(w/pi,abs(H123(1:length(w)))/max(abs(H123(1:length(w))))*max(abs(XNOISE(1:length(w)))),'r');
xlim([0 1]);
title('Espectro de la se�al original');
xlabel('\omega (rad/s)');
ylabel('Amplitud');
str1 = 'Se�al original';
str2 = 'Respuesta en m�dulo del filtro '; 
legend(str1,str2,'Location','West');

subplot(2,1,2);
plot(w/pi,abs(Y(1:length(w))),'r');
xlim([0 1]);
grid on;
title('Espectro de la se�al filtrada')
xlabel('\omega (rad/s)')
ylabel('Amplitud')

% audio = wavread(pista_01.wav);
% audio_noise = wavread(pista_01_noise.wav);
% audio_cleaned = wavread(pista_01_cleaned.wav);

%sound("pista_01.wav",fs);
%sound(pista_01_noise,fs);
%sound(pista_01_cleaned,fs);


% figure();
% plot(w,20*log10(abs(H(1:length(w)))-1))

% grid on;
% title('H(jw) - 1, para mirar la  amplitud del ripple');

% Se exporta el filtro notch fir

save notch_fir.mat h123;
