close all,clc

% Cargar la imagen
imagen = imread('imagen2.jpg');

% Obtener las componentes de color R, G y B
R = imagen(:, :, 1);
G = imagen(:, :, 2);
B = imagen(:, :, 3);

%Normalizar las componentes de color (opcional)
 R = double(R) / 255;
 G = double(G) / 255;
 B = double(B) / 255;

% Convertir las componentes de color en secuencias binarias
umbral = 0.5; % Puedes ajustar este umbral según tus necesidades
secuencia_binaria_R = R > umbral;
secuencia_binaria_G = G > umbral;
secuencia_binaria_B = B > umbral;

% Combinar las secuencias binarias en una sola secuencia
secuencia_total = [secuencia_binaria_R(:); secuencia_binaria_G(:); secuencia_binaria_B(:)];

% Transmitir la secuencia binaria como una secuencia de dígitos binarios
% Puedes enviar la secuencia total directamente al receptor.
secuencia_binaria_transmitida = reshape(secuencia_total, 1, []);


% Parámetro M para la modulación M'-QAM
M = 256; % Puedes cambiar M según tus necesidades

% Genera una constelación M'-QAM de forma automática
constelacion = qammod(0:M-1, M, 'gray');

% Número de bits por símbolo en la modulación M'-QAM
bits_por_simbolo = log2(M);

% Secuencia de bits de entrada (ejemplo)
bits_recibidos = secuencia_binaria_transmitida;
tamano_bits=length(bits_recibidos);
%Calcula si se deben agregar 0 al vector si no cumple con las dimensiones
%correctas
modulo = mod(tamano_bits,bits_por_simbolo);
disp(['Modulo:',num2str(modulo)]);
% La longitud de bits debe ser un múltiplo de bits_por_simbolo
ceros_final=0;
if modulo~= 0
    ceros_final=bits_por_simbolo-modulo;
    disp(['Ceros agregados:',num2str(ceros_final)]);    
    % Crear un vector de ceros del tamaño deseado
    cerosAgregados = zeros(1, ceros_final);

    % Concatenar el vector de ceros con tu vector original
    bits = [bits_recibidos, cerosAgregados];
    % Eliminar los últimos n elementos para obtener el vector original:Usar
    % en Recepcion
    %vectorOriginal = bits(1:end-ceros_final);
else
    bits=bits_recibidos;
end

% Divide los bits en grupos de bits_por_simbolo
grupos_de_bits = reshape(bits, bits_por_simbolo, length(bits)/bits_por_simbolo).';

% Modula los grupos de bits utilizando la constelación
msbfirst = true;
valores_alfabeto = constelacion(bi2de(grupos_de_bits, 'left-msb') + 1);

% Imprime la constelación generada
% disp('Simbolos generados:');
% disp(valores_alfabeto);

% Obtiene partes reales e imaginarias utilizando las funciones real e imag
partes_reales = real(valores_alfabeto);
partes_imaginarias = imag(valores_alfabeto);

% Número de ceros a insertar entre cada elemento
n = 8;

% Inicializa un vector para almacenar el resultado
%Anade n+1 ceros en la funcion original para convolucionarla
senal_modulada_imag = upsample(partes_imaginarias,n+1);
senal_modulada_real = upsample(partes_reales,n+1);

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

Es = sum((abs(constelacion).^2)./M);
% Eb=Es/(log2(M));
EbNo=25;
ebno = 10^((EbNo)/10);
sigma = sqrt(Es/(2*log2(M)*ebno));


Z = sigma.*randn(1,length(X));%Ruido AWGN

%Demodulacion
Y=X+Z; %Senal recibida despues del canal AWGN
%Baja la senal de pasa banda a banda base
y_real=sqrt(2)*Y.*cos(2*pi*fc.*t);
y_imag=-sqrt(2)*Y.*sin(2*pi*fc.*t);


%Quitar los ceros a la senal

%Filtra la senal resultante con el filto p(-t)
y_prima_real=filter(h,1,y_real);
y_prima_imag=filter(h,1,y_imag);

%Toma muestras de la senal cada n+1 muestras con un desplazamiento de fase n 
y_prima_real_sample = downsample(y_prima_real,n+1,n);

y_prima_imag_sample = downsample(y_prima_imag,n+1,n);

%Unir la parte real e imaginaria recibida para recuperar nuevamente los
%simbolos
simbolos_recibidos=complex(y_prima_real_sample,y_prima_imag_sample);
% disp('Simbolos recibidos')
% disp(simbolos_recibidos);


% Inicializar variables para almacenar la distancia mínima y los índices correspondientes
distanciaMinima = Inf;
indicesMinimos = [];
for j = 1:length(simbolos_recibidos)
    for i = 1:length(constelacion)
        distancia = abs(simbolos_recibidos(j) - constelacion(i));
        
        % Verificar si la distancia actual es menor que la distancia mínima
            if distancia < distanciaMinima
                    distanciaMinima = distancia;
                    indicesMinimos(j) = i;
            end
    end
            distanciaMinima = Inf;
end


% Mostrar los resultados
vector_distancias=[];
for j = 1:length(simbolos_recibidos)
    vector_distancias(j)=constelacion(indicesMinimos(j));
end
%Asigna un simbolo de la constelacion a los simbolos recibidos (decision)
simbolos_decision=vector_distancias;

%Realiza el demapeo de los simbolos recibidos en Bits
demap= qamdemod(simbolos_decision,M,'gray');
bits_bin=dec2bin(demap);
bits_matrix=reshape(bits_bin.', 1, []);
bits_demap_transpose=str2num(bits_matrix.');
bits_demap=transpose(bits_demap_transpose);

if ceros_final~= 0
    posicion_ceros = length(bits_demap);
    posicion_eliminar = max(1, posicion_ceros-ceros_final);
    % Trunca el vector hasta la posición de eliminación
    vector_sin_ultimos_ceros = bits_demap(1:posicion_eliminar);
    secuencia_binaria_recibida=vector_sin_ultimos_ceros;
else
    secuencia_binaria_recibida=bits_demap;
end

%Compara la secuencia recibida y la transmitida
vector_tx=double(secuencia_binaria_transmitida);
comparacion=symerr(secuencia_binaria_recibida,vector_tx)

% Recuperacion de la imagen en Recepcion

% Paso 1: Separar la secuencia binaria en las tres componentes de color
num_elementos = numel(secuencia_binaria_recibida) / 3; % Número de elementos por componente
secuencia_rec_binaria_R = secuencia_binaria_recibida(1:num_elementos);
secuencia_rec_binaria_G = secuencia_binaria_recibida(num_elementos+1:2*num_elementos);
secuencia_rec_binaria_B = secuencia_binaria_recibida(2*num_elementos+1:end);

% Paso 2: Reconstruir las componentes de color (inversa de la binarización)
umbral_inverso = 0.5; % Este es el umbral inverso
imagen_reconstruida_R = reshape(secuencia_rec_binaria_R, size(imagen, 1), size(imagen, 2)) * 255;
imagen_reconstruida_R = uint8(imagen_reconstruida_R);

imagen_reconstruida_G = reshape(secuencia_rec_binaria_G, size(imagen, 1), size(imagen, 2)) * 255;
imagen_reconstruida_G = uint8(imagen_reconstruida_G);

imagen_reconstruida_B = reshape(secuencia_rec_binaria_B, size(imagen, 1), size(imagen, 2)) * 255;
imagen_reconstruida_B = uint8(imagen_reconstruida_B);

% Paso 3: Combina las tres componentes de color (R, G y B) en una imagen a color
imagen_reconstruida = cat(3, imagen_reconstruida_R, imagen_reconstruida_G, imagen_reconstruida_B);


% %Grafica la constelación generada
figure(1);
scatter(real(valores_alfabeto), imag(valores_alfabeto), 'filled');
title(['Diagrama de Constelación ',num2str(M), '-QAM']);
xlabel('Parte Real');
ylabel('Parte Imaginaria');
grid on;
%axis([-4 4 -4 4]); % Ajusta los límites del gráfico según sea necesario

%Grafica de las senales Modulada y Demodulada
% figure(2);
% subplot(2,2,1);
% stem(conformed_pulses_real)
% xlim([0,180]);
% title('Parte Real multiplicada por p(t)')
% 
% subplot(2,2,2);
% stem(conformed_pulses_imag)
% xlim([0,180]);
% title('Parte Imaginaria multiplicada por p(t)')
% subplot(2,2,3)
% stem(y_prima_real_sample)
% title('Simbolos Reales Recibidos')
% xlim([0,20]);
% subplot(2,2,4)
% stem(y_prima_imag_sample)
% title('Simbolos Imaginarios Recibidos')
% xlim([0,20]);
%Grafico de los simbolos recuperados.

figure(3);
scatter(y_prima_real_sample, y_prima_imag_sample, 'filled');
title(['Diagrama de Constelación Recibida',num2str(M), '-QAM']);
xlabel('Parte Real');
ylabel('Parte Imaginaria');
grid on;
%axis([-4 4 -4 4]); % Ajusta los límites del gráfico según sea necesario

figure(4)
subplot(2,1,2);
imshow(imagen_reconstruida);
title('Imagen Reconstruida')
subplot(2,1,1);
imshow(imagen);
title('Imagen Original')

figure(5)
scatter(real(simbolos_recibidos), imag(simbolos_recibidos), 'filled');
hold on
scatter(real(constelacion), imag(constelacion), 'filled');
legend('Simbolos Recibidos','Constelacion Original')
grid on
hold off