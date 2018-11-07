%=======================================================================
% Procesamiento de Señales en Tiempo Discreto
% Prof.: Dr. Mario Hueda
% Practico Lab. 2
%=======================================================================
clear;
close all;

%============================================
% Generacion de la Respuesta al Impulso
%============================================
fB = 32e9;	% Velocidad de simbolos (baud rate)
T = 1/fB; % Tiempo entre simbolos
M = 8;  %Factor de sobremuestreo
fs = fB*M;	% Sample rate

beta = 0.5001; %Factor de roll-off
L = 10;  % 2*L*M+1 es el largo del filtro sobremuestreado
t = [-L:1/M:L]*T;
n_delay_RC_filter = L*M; %Retardo del filtro RC
gn = sinc(t/T).*cos(pi*beta*t/T)./(1-4*beta^2*t.^2/T^2);

Lf=100;
n=[-Lf:Lf];
n_delay_Hilbert_filter = Lf; %Retardo del filtro de Hilbert
fn = 2*sin(pi*n/2).^2./(pi*n);
fn(Lf+1)=0;

figure(1)
subplot 211
h = stem(gn);
title('Respuesta al Impulso del Filtro RC');
xlabel('n');
grid on
subplot 212
h = stem(fn);
title('Respuesta al Impulso del Transformador de Hilbert');
xlabel('n');
grid on

%break

%============================================
% Calculo de la Respuesta en Frecuencia
%============================================
Omega = [-1:1/2^8:1]*pi;
N = 1000;
n = [0:N];
index = 1;
for omega=Omega
    xn = exp(j*omega*n);
    yn = conv(xn,fn);
    H_Mag(index) = abs(yn(N/2));
    H_Fase(index) = angle(yn(N/2)*conj(xn(N/2-n_delay_Hilbert_filter)));
    index = index+1;
end

figure(2)
subplot 211
h=plot(Omega/pi,H_Mag);
title('Magnitud');
ylabel('|H|')
xlabel('\Omega/\pi');
grid
subplot 212
h=plot(Omega/pi,H_Fase);
title('Fase');
ylabel('angle(H)')
xlabel('\Omega/\pi');
grid

%break

%============================================
% Generacion Simbolos
%============================================
n_symbols = 2^14;
ak = 2*randint(1,n_symbols)-1+j*(2*randint(1,n_symbols)-1);
xn = zeros(1,n_symbols*M);
xn(1:M:end) = ak;

figure(3)
h = plot(ak,'x');
xlabel('Real(ak)')
ylabel('Imag(ak)');
axis([-1.2 1.2 -1.2 1.2])
%break

%============================================
% Señal Banda-Base
%============================================

yn = conv(xn,gn);
figure(4)
q = spectrum.welch;      
Hpsd = psd(q,yn,'nfft',1024);
h=plot(Hpsd);
title('PSD Señal en Banda Base')

%============================================
% Señal Modulada (Analitica)
%============================================

n=[1:length(yn)];
Omega_c=3*pi/M; %Portadora arbitraria para generar señal analitica (ojo: depende de M)
carrier=exp(j*Omega_c*n);
sn = yn.*carrier;

figure(5)
q = spectrum.welch; 
Hpsd = psd(q,sn,'nfft',1024);
h=plot(Hpsd);
title('PSD Señal Modulada (Analitica)')

%============================================
% Señal Transmitida (Parte Real de la Señal Analitica)
%============================================

sn_r = real(sn);
figure(6)
q = spectrum.welch;
Hpsd = psd(q,sn_r,'nfft',1024,'SpectrumType','twosided');
h=plot(Hpsd);
title('PSD Señal Transmitida (Real)')

%============================================
% Filtro de Particion de Fase
%============================================

sn_i = conv(sn_r,fn);
sn_hat = sn_r+j*sn_i(Lf+1+0:end-Lf+0);
figure(7)
q = spectrum.welch;
Hpsd = psd(q,sn_hat,'nfft',1024,'SpectrumType','twosided');
h=plot(Hpsd);
title('PSD Señal a la Salida del Filtro de Particion de Fase')

%============================================
% Comparacion Partes Imaginarias (transmitida y recibida) 
%============================================

sn_i = conv(sn_r,fn);
sn_hat = sn_r+j*sn_i(Lf+1:end-Lf);
figure(8)
n=[1:100]+1000; %100 puntos de ventana de tiempo arbitraria
h=plot(n,imag(sn(n)),'b',n,sn_i(n+Lf),'ro');
legend('Transmitida', 'Recuperada');
title('Parte Imaginaria de la Señal Analitica')

%============================================
% Señal Demodulada (Banda-Base)
%============================================

n=[1:length(sn_hat)];
carrier=exp(j*Omega_c*n);
yn_hat = sn_hat.*conj(carrier); %portadora conjugada (demodulacion)
figure(9)
q = spectrum.welch;
Hpsd = psd(q,yn_hat,'nfft',1024,'SpectrumType','twosided');
h=plot(Hpsd);
title('PSD Señal Demodulada (Banda-Base)')

