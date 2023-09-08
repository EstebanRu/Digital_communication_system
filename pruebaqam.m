clear all, close all
% Parámetro M para la modulación M'-QAM
M = 4; % Puedes cambiar M según tus necesidades

% Genera una constelación M'-QAM de forma automática
constelacion = qammod(0:M-1, M, 'gray');

% Número de bits por símbolo en la modulación M'-QAM
bits_por_simbolo = log2(M);

% Secuencia de bits de entrada (ejemplo)
bits = [0 0 0 1 1 0 1 1 1 0];

% La longitud de bits debe ser un múltiplo de bits_por_simbolo
if mod(length(bits), bits_por_simbolo) ~= 0
    error('La longitud de la secuencia de bits no es un múltiplo de %s.',int2str(bits_por_simbolo));
end

% Divide los bits en grupos de bits_por_simbolo
grupos_de_bits = reshape(bits, bits_por_simbolo, length(bits)/bits_por_simbolo).';

% Modula los grupos de bits utilizando la constelación
msbfirst = true;
valores_alfabeto = constelacion(bi2de(grupos_de_bits, 'left-msb') + 1);

% Imprime la constelación generada
disp('Constelación generada:');
disp(valores_alfabeto);

% Obtiene partes reales e imaginarias utilizando las funciones real e imag
partes_reales = real(valores_alfabeto);
partes_imaginarias = imag(valores_alfabeto);

% Número de ceros a insertar entre cada elemento
n = 8;

% Inicializa un vector para almacenar el resultado
senal_modulada_imag = [];

% Itera sobre el vector original y agrega los ceros entre cada elemento
for i = 1:length(partes_imaginarias)
    senal_modulada_imag = [senal_modulada_imag, partes_imaginarias(i), zeros(1, n)];
end

senal_modulada_real = [];
for i = 1:length(partes_reales)
    senal_modulada_real = [senal_modulada_real, partes_reales(i), zeros(1, n)];
end

%Implementa el filtro conformador p(t) PARTE REAL
h = rcosdesign(0.5, 1, n, 'sqrt');
conformed_pulses_real= filter(h,1,senal_modulada_real);

%Implementa el filtro conformador p(t) PARTE IMAGINARIA
conformed_pulses_imag= filter(h,1,senal_modulada_imag);

%Suma la Parte Real e Imaginaria para salir por el canal
senal_suma= conformed_pulses_imag+conformed_pulses_real;


%Grafica la constelación
figure(1);
scatter(real(valores_alfabeto), imag(valores_alfabeto), 'filled');
title(['Diagrama de Constelación ',num2str(M), '-QAM']);
xlabel('Parte Real');
ylabel('Parte Imaginaria');
grid on;
axis([-4 4 -4 4]); % Ajusta los límites del gráfico según sea necesario

%Muestra la secuencia de simbolos convolcuinada con el filtro conformador
figure(2);
subplot(2,1,1);
stem(conformed_pulses_imag)
title('Parte Imaginaria Modulada')
subplot(2,1,2);
stem(conformed_pulses_real)
title('Parte Real Modulada')

