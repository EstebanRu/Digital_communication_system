%Cargar la imagen
clear all, clc , close all;
ruta = 'C:\Users\Pablo Garzon\Desktop\Matlab\coin.jpg';
imagen = imread(ruta);
%convertir imagen en ua secuencia binaria
datos_binarios = reshape(de2bi(imagen), 1 , []);
% Datos binarios de entrada (asegúrate de que la longitud sea un múltiplo de 4)


% Definir el mapeo de símbolos 16-QAM (I y Q)
simbolos_I = [-3 -3 -3 -3 -1 -1 -1 -1 3 3 3 3 1 1 1 1];
disp(simbolos_I(1))
simbolos_Q = [-3 -1 1 3 -3 -1 1 3 -3 -1 1 3 -3 -1 1 3];

% Inicializar el vector de símbolos 16-QAM
simbolos_16QAM = zeros(1, length(datos_binarios)/4);
% Mapear grupos de 4 bits a símbolos 16-QAM
for i = 1:length(simbolos_16QAM)
    indice = bin2dec(num2str(datos_binarios((i-1)*4 + 1:i*4))) + 1; % +1 para convertir a índice
    simbolos_16QAM(i) = simbolos_I(indice) + 1i * simbolos_Q(indice);
end

% Visualizar los símbolos 16-QAM resultantes
scatterplot(simbolos_16QAM);
title('Diagrama de Constelación 16-QAM');

% Parámetros
SNR_dB = 20;                  % Relación señal-ruido (SNR) en dB
num_antenas_tx = 1;           % Número de antenas en el transmisor
num_antenas_rx = 1;           % Número de antenas en el receptor

% Simulación del canal Rayleigh (multipath)
canal_rayleigh = (randn(num_antenas_rx, num_antenas_tx) + 1i * randn(num_antenas_rx, num_antenas_tx)) / sqrt(2);
simbolos_rx_multipath = canal_rayleigh * simbolos_16QAM;

% Simulación del ruido AWGN
SNR = 10^(SNR_dB/10);  % Convertir SNR de dB a escala lineal
potencia_senal = mean(abs(simbolos_rx_multipath).^2);
potencia_ruido = potencia_senal / SNR;
ruido = sqrt(potencia_ruido) * (randn(size(simbolos_rx_multipath)) + 1i * randn(size(simbolos_rx_multipath)));

% Señal recibida en el receptor (con ruido)
simbolos_rx_con_ruido = simbolos_rx_multipath + ruido;

% Demodulación de los símbolos recibidos
simbolos_demodulados = qamdemod(simbolos_rx_con_ruido, 16);

% Calcular la tasa de error de bits (BER)
bits_tx = de2bi(qamdemod(simbolos_16QAM, 16), 4);
bits_rx = de2bi(simbolos_demodulados, 4);
errores = sum(sum(bits_tx ~= bits_rx));
BER = errores / ( length(datos_binarios)/4 * log2(16));

fprintf('Tasa de Error de Bits (BER): %f\n', BER);

