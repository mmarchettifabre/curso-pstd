%=======================================================================
% Procesamiento de Señales en Tiempo Discreto
% Prof.: Dr. Mario Hueda
% Practico Lab. 5
%=======================================================================
clear all;
close all;
ON=1;
OFF=0;
f0=1.0e6;              % Ancho de banda, en Hz
M=2^4;                 % Factor de sobremuestreo del SD (sigma-delta)
R=4;                   % Factor de sobremuestreo despues del diezmado (>2 para cumplir Nyquist)
fs=R*f0*M;             % Frecuencia de muestreo, in Hz
Ts=1/fs;
N=15000;


%============================================
% Generacion de senal en bandabase
%============================================
g=rcosine(f0,fs,'normal',.1,40);  % Filtro Tx
ak = 2*(randi(2,1,N)-1)-1;
xn = zeros(1,N*fs/f0);
xn(1:fs/f0:end) = ak;
x=.125*filter(g,1,xn);

%============================================
% Conversor Sigma-Delta de 2do Orden y Filtro CIC de 2 Etapas
%============================================
e1=zeros(1,N*M);
e2=zeros(1,N*M);
s1=zeros(1,N*M);
s2=zeros(1,N*M);
y=ones(1,N*M);
xhat=zeros(1,N);
xhat1=zeros(1,N);
xhat2=zeros(1,N);
k=2;
y_acum1=0;
y_acum2=0;
for n=2:N*M
    e1(n)=x(n)-y(n-1);
    s1(n)=e1(n)+s1(n-1);
    e2(n)=s1(n)-y(n-1);
    s2(n)=e2(n)+s2(n-1);
    y(n)=(s2(n)>0)-(s2(n)<=0);
    %y(n)=s2(n)*1.0+4.*(rand-.5); % Prueba con ruido blanco
    y_acum1=mod(y_acum1+y(n)+1,2*M^2); % La senal de entrada del CIC es positiva [0 2] con media 1!
    y_acum2=mod(y_acum2+y_acum1,2*M^2);
    if mod(n,M)==0
        xhat1(k)=y_acum2;
        xhat2(k)=mod(xhat1(k)-xhat1(k-1),2*M^2);
        xhat(k)=mod(xhat2(k)-xhat2(k-1),2*M^2);
        k=k+1;
    end
end
xhat=xhat(2000:end)/M^2-1; % Ajuste de ganancia (DC=1) y eliminacion de continua
figure(1)
subplot 211
[H F]=psd(y(2000:end),2^10,fs);
h=plot(F,10*log10(H),F,20*log10(abs((1-exp(-j*2*pi*F*Ts)).^2)));
set(h,'Linewidth',2);
set(h,'Markersize',16);
set(gca,'XScale','lin','YScale','lin','FontWeight','bold','FontSize',14);
set(gca,'Linewidth',2);
xlabel('f [Hz]');
ylabel('PSD');
grid on

subplot 212
[H F]=psd(xhat,2^10,fs/M);
h=plot(F,10*log10(H/max(H)));
set(h,'Linewidth',2);
set(h,'Markersize',16);
set(gca,'XScale','lin','YScale','lin','FontWeight','bold','FontSize',14);
set(gca,'Linewidth',2);
xlabel('f [Hz]');
ylabel('PSD');
grid on
%break

%============================================
% Diagramas Ojos de la Senal Diezmada
%============================================
eyediagram(xhat,R,1,2)

%============================================
% Diagramas Ojos de la Senal Diezmada y Filtrada
%============================================
n=[-30:1.25/R:30];
f=sinc(n); % Filtro pasabajo (saca ruido de cuantizacion fuera de la banda de la senal)
xhat_f=filter(f/sum(f),1,xhat);
eyediagram(xhat_f(50*R:end),R,1,3)

%============================================
% PSD de la Senal Diezmada antes y despues del LPF
%============================================
figure
[H F]=psd(xhat,2^10,fs/M);
[Hf F]=psd(xhat_f,2^10,fs/M);
h=plot(F,10*log10(H/max(H)),F,10*log10(Hf/max(Hf)));
set(h,'Linewidth',2);
set(h,'Markersize',16);
set(gca,'XScale','lin','YScale','lin','FontWeight','bold','FontSize',14);
set(gca,'Linewidth',2);
xlabel('f [Hz]');
ylabel('PSD');
grid on
