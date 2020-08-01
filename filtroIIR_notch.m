clear all;
close all;

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

noise = A1*cos(2*pi*f1*t) + A2*cos(2*pi*f2*t) + A3*cos(2*pi*f3*t) ;
[h,w]=freqz(noise,8192);
% figure, plot(w,abs(h),'b');
% 
% figure,plot(t,noise,'r');
% xlim([0 0.08]);
% grid on;
% str='0.08s de interferencia';
% title(str);
% auxFFT  = abs(fft(noise,2^nextpow2(length(noise))));
% w       = 0:2*pi/(length(auxFFT)-1):2*pi;
% figure,plot(w(1,1:floor(length(w)/2)),auxFFT(1,1:floor(length(w)/2)),'r');
% grid on;
% str='FFT de la interferencia';
% title(str);

xconruido = x + noise';

%% Filtro IIR butterworth
% Realizo primero un filtro pasa bajos
% Primero calculo D y K en funcion de ds, dp,wp y ws
% Parametros
dp=0.1;
ds=0.06;
wp=0.016*pi; %frec inferior
ws=0.020*pi; %frec superior
T=1;
wap = wp/T;
was = ws/T;
Ap = -20*log10(1-dp);
As = -20*log10(ds);

K = wp/ws;
%% Calculo de N y wc
N = 4
wc = 0.017*pi;
[z,p,k] = buttap(N); 
num = real(poly(z)); % Numerador de la transferencia
num = num*(wc^N)*k;
den= real(poly(p*wc)); % Denominador de la transferencia
[numd,dend] = impinvar(num,den,T);
sys =tf(numd,dend);
x=filter(numd,dend,x); % Filtro

w =(1:length(x)/2)*2*pi/length(x);
H = freqz(numd,dend,w); 
Hma = abs(H);
Hfase = angle(H);
Hmdb = 20*log10((Hma+eps)/(max(Hma)));
% figure,plot(w, Hma);
% figure,plot(w, Hfase)

[numsb,densb] = iirlp2bs(numd,dend,0.0534,[0.0503 0.0628]);
h1 =tf (numsb,densb);
Hsb = freqz(numsb,densb,w); 
Hmasb = abs(Hsb);
Hfasesb = angle(Hsb);
Hmdbsb = 20*log10((Hmasb+eps)/(max(Hmasb)));
% figure,plot(w, Hmasb);
% figure,plot(w, Hfasesb);

%fvt = fvtool(numd,dend,num,den);
%legend(fvt,'Prototype','Target')

%% Filtro IIR butterworth
% Realizo primero un filtro pasa bajos
% Primero calculo D y K en funcion de ds, dp,wp y ws
% Parametros
dp=0.1;
ds=0.06;
wp2=0.0090*pi; %frec inferior
ws2=0.01*pi; %frec superior
T=1;
wap = wp2/T;
was = ws2/T;
Ap = -20*log10(1-dp);
As = -20*log10(ds);

K = wp/ws;
%% Calculo de N y wc
N = 3
wc2 = 0.0095*pi;
[z2,p2,k2] = buttap(N); 
num2 = real(poly(z2)); % Numerador de la transferencia
num2 = num2*(wc2^N)*k2;
den2= real(poly(p2*wc2)); % Denominador de la transferencia
[numd2,dend2] = impinvar(num2,den2,T);
sys =tf(numd2,dend2);
x=filter(numd2,dend2,x); % Filtro

w =(1:length(x)/2)*2*pi/length(x);

H2 = freqz(numd2,dend2,w); 
Hma2 = abs(H2);
Hfase2 = angle(H2);
Hmdb2 = 20*log10((Hma2+eps)/(max(Hma2)));
% figure,plot(w, Hma2);
% figure,plot(w, Hfase2);

[numsb2,densb2] = iirlp2bs(numd2,dend2,0.0298,[0.0258 0.0328]);
h2 =tf (numsb2,densb2);
Hsb2 = freqz(numsb2,densb2,w); 
Hmasb2 = abs(Hsb2);
Hfasesb2 = angle(Hsb2);
Hmdbsb2 = 20*log10((Hmasb2+eps)/(max(Hmasb2)));
% figure,plot(w, Hmasb2);
% figure,plot(w, Hfasesb2);

%% Filtro IIR butterworth
% Realizo primero un filtro pasa bajos
% Primero calculo D y K en funcion de ds, dp,wp y ws
% Parametros
dp=0.1;
ds=0.06;
wp=0.0305*pi; %frec inferior
ws=0.0325*pi; %frec superior
T=1;
wap = wp/T;
was = ws/T;
Ap = -20*log10(1-dp);
As = -20*log10(ds);

K = wp/ws;
%% Calculo de N y wc
N = 4
wc3 = 0.1025;
[z3,p3,k3] = buttap(N); 
num3 = real(poly(z3)); % Numerador de la transferencia
num3 = num3*(wc3^N)*k3;
den3= real(poly(p*wc3)); % Denominador de la transferencia
[numd3,dend3] = impinvar(num3,den3,T);
sys =tf(numd3,dend3);
x=filter(numd3,dend3,x); % Filtro

w =(1:length(x)/2)*2*pi/length(x);
H3 = freqz(numd3,dend3,w); 
Hma3 = abs(H3);
Hfase3 = angle(H3);
Hmdb3 = 20*log10((Hma3+eps)/(max(Hma3)));
%figure,plot(w, Hma3);
%figure,plot(w, Hfase3);

[numsb3,densb3] = iirlp2bs(numd3,dend3,0.1025,[0.1 0.1067]);
h3 = tf (numsb3,densb3);
h3d = c2d(h3,1);
Hsb3 = freqz(numsb3,densb3,w); 
Hmasb3 = abs(Hsb3);
Hfasesb3 = angle(Hsb3);
Hmdbsb3 = 20*log10((Hmasb3+eps)/(max(Hmasb3)));
%figure,plot(w, Hmasb3);
%figure,plot(w, Hfasesb3);

h1f = dfilt.df2(numsb,densb);
h2f = dfilt.df2(numsb2,densb2);
h3f = dfilt.df2(numsb3,densb3);

hd = dfilt.cascade(h1f, h2f,h3f);
HD = freqz(hd,w);

HDnotch = abs(HD);
HDfase = angle(HD);
figure,plot(w, HDnotch);
figure,plot(w, HDfase);

%hfvt= fvtool(hd,'Color','white');
