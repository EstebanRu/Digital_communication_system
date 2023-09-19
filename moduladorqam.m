clear all, close all
% Parámetro M para la modulación M'-QAM
M = 16; % Puedes cambiar M según tus necesidades

% Genera una constelación M'-QAM de forma automática
constelacion = qammod(0:M-1, M, 'gray');

% Número de bits por símbolo en la modulación M'-QAM
bits_por_simbolo = log2(M);

% Secuencia de bits de entrada (ejemplo)
bits = [0 0 1 1 1 1 1 0 0 0 1 0 0 1 0 1 0 0 0 1 1 0 1 1 1 0 1 0];
tamano_bits=length(bits)

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

%Configuracion de velocidad de transmision, factor y frecuencia de muestreo
Rs=100;
U=6;
fs=U*Rs;
ts=1/fs;
fc=2*Rs;
Ns=length(conformed_pulses_real);

%Multiplicacion por las portadoras
t= (0:Ns-1)/fs;
x_real=sqrt(2)*conformed_pulses_real.*cos(2*pi*fc.*t);
x_imag=-sqrt(2)*conformed_pulses_imag.*sin(2*pi*fc.*t);
X=x_imag + x_real;

%Demodulacion
Y=X;
%Baja la senal de pasa banda a banda base
y_real=sqrt(2)*Y.*cos(2*pi*fc.*t);
y_imag=-sqrt(2)*Y.*sin(2*pi*fc.*t);


%Quitar los ceros a la senal

%Filtra la senal resultante con el filto p(-t)
y_prima_real=filter(h,1,y_real);
y_prima_imag=filter(h,1,y_imag);

%Toma muestras de la senal cada n muestras
y_prima_real_sample = [];

% Usar un bucle for para recorrer el vector original y seleccionar cada quinta muestra
for i = n+1:n+1:length(y_prima_real)
    muestra = y_prima_real(i);
    y_prima_real_sample = [y_prima_real_sample, muestra];
end

%Toma muestras de la senal cada n muestras
y_prima_imag_sample = [];

% Usar un bucle for para recorrer el vector original y seleccionar cada quinta muestra
for i = 9:n+1:length(y_prima_imag)
    muestra_imag = y_prima_imag(i);
    y_prima_imag_sample = [y_prima_imag_sample, muestra_imag];
end
figure(4)
subplot(2,1,1)
stem(y_prima_real_sample)
title('Simbolos Reales')
subplot(2,1,2)
stem(y_prima_imag_sample)
title('Simbolos Imaginarios')

%Grafica la constelación
figure(1);
scatter(real(valores_alfabeto), imag(valores_alfabeto), 'filled');
title(['Diagrama de Constelación ',num2str(M), '-QAM']);
xlabel('Parte Real');
ylabel('Parte Imaginaria');
grid on;
axis([-4 4 -4 4]); % Ajusta los límites del gráfico según sea necesario

%Grafica de las senales Modulada y Demodulada
figure(2);
subplot(2,2,1);
stem(conformed_pulses_real)
title('Parte Real multiplicada por p(t)')
subplot(2,2,2);
stem(conformed_pulses_imag)
title('Parte Imaginaria multiplicada por p(t)')

subplot(2,2,3)
stem(y_prima_real);
title('Senal despues de p(-t) Real')
subplot(2,2,4)
stem(y_prima_imag);
title('Senal despues de p(-t) Imaginaria')

%Grafico de la senal multiplicada por la portadora
figure(3);
subplot(2,1,1)
stem(Y);
title('Senal Transmitida')
subplot(2,2,3);
stem(y_real)
title('Parte Real de y(t)')
subplot(2,2,4);
stem(y_imag)
title('Parte Imaginaria de y(t)')


