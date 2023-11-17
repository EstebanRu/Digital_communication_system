close all,clear all,clc

nombre_archivo = 'imagen2.jpg'; % Ruta de la imagen
imagen_rgb = imread(nombre_archivo); % Cargar la imagen RGB

umbral=128;

% Extraer los canales de color rojo, verde y azul
canal_rojo = imagen_rgb(:,:,1);
canal_verde = imagen_rgb(:,:,2);
canal_azul = imagen_rgb(:,:,3);
[filas0,columnas0]=size(canal_rojo);

vector_rojo = reshape(canal_rojo, 1, []);
vector_verde = reshape(canal_verde, 1, []);
vector_azul = reshape(canal_azul, 1, []);

vector_r = vector_rojo > umbral;
vector_v = vector_verde > umbral;
vector_a = vector_azul >umbral;
[filas,columnas]=size(vector_r);

vectorCompleto = horzcat(vector_r, vector_v,vector_a);

secuencia_binaria_transmitida = vectorCompleto;

% Parámetro M para la modulación M'-QAM
M = 16; % Puedes cambiar M según tus necesidades

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
simbolos_generados = constelacion(bi2de(grupos_de_bits, 'left-msb') + 1);

dura_sim=0:0.1:length(simbolos_generados);
% Obtiene partes reales e imaginarias utilizando las funciones real e imag
partes_reales = real(simbolos_generados);
partes_imaginarias = imag(simbolos_generados);
%Configuracion de velocidad de transmision, factor y frecuencia de muestreo
Rs=100;
U=6;
fs=U*Rs;
ts=1/fs;
fc=2*Rs;
%% 
beta=0.49999;
T=1;
S=8;
tp=-S:1/S:S;
alpha=[1 0.9 0.2];
tau=[0 T 2*T];

r=zeros(1,length(tp));

for i=1:length(alpha)
   r=r+alpha(i)*exp(-1j*2*pi*fc*tau(i))*(sinc((tp-tau(i))/T).*cos((pi*beta*(tp-tau(i)))/T))./(1-((2*beta*(tp-tau(i)))/T).^2);
end
figure,plot(tp,real(r),'b',tp,imag(r),'r'),grid,xlabel('t/T'),ylabel('r(t)'),title('Respuesta global al impulso')

rk=downsample(r,S);
figure();
stem(real(rk),'.b'),hold on,stem(imag(rk),'.r'),grid,xlabel('t/T'),ylabel('r(kT)'),title('Canal discreto equivalente')
%Convolucion simbolos con el canal discreto equivalente
uk=conv(simbolos_generados,rk);
figure();
subplot(211),stem(real(simbolos_generados),'.r'),hold on,stem(imag(simbolos_generados),'.b'),title('Simbolos Transmitidos'),legend('simbolos reales','simbolos imag'),xlim([0,20]),grid on;
subplot(212),stem(real(uk),'.r'),hold on,stem(imag(uk),'.b'),grid,title('uk'),legend('simbolos reales','simbolos imag'),xlim([S+1,S+20]);
figure();
scatter(real(simbolos_generados),imag(simbolos_generados),'filled'),grid,axis('equal'),xlabel('Real'),ylabel('Imag'),title(['Constelacion ',num2str(M), ' -QAM']);
figure();
scatter(real(uk(S+1:S+length(simbolos_generados))),imag(uk(S+1:S+length(simbolos_generados))),'filled'),grid,axis('equal'),xlabel('Real'),ylabel('Imag'),title('Constelacion recibida con multitrayecto')

 
%% CANAL DISCRETO EQUIVALENTE + RUIDO
%RUIDO COMPLEJO AWGN
Es = sum((abs(constelacion).^2)./M);
% Eb=Es/(log2(M));
EbNo=20;
ebno = 10^((EbNo)/10);
sigma = sqrt(Es/(2*log2(M)*ebno));
Z = sigma.*(randn(1,length(uk))+j*randn(1,length(uk)));%Ruido AWGN

uk=uk+Z; %Senal recibida despues del canal AWGN

figure();
stem(real(uk(S+1:S+length(simbolos_generados))),'.r'),hold on,stem(imag(uk(S+1:S+length(simbolos_generados))),'.b'),grid,title('Observaciones uk, sin ecualizacion'),legend('simbolos reales','simbolos imag'),xlim([0,S+20]);
figure();
scatter(real(simbolos_generados),imag(simbolos_generados),'filled'),grid,axis('equal'),xlabel('Real'),ylabel('Imag'),title('Constelacion generada')
figure();
scatter(real(uk(S+1:S+length(simbolos_generados))),imag(uk(S+1:S+length(simbolos_generados))),'filled'),grid,axis('equal'),xlabel('Real'),ylabel('Imag'),title('Constelacion recibida con ruido y multitrayecto')

%% ECUALIZADOR ZF
 valoresR = rk(S+1:S+length(alpha));

figure()
subplot(211),stem(real(uk)),title('Observaciones sin ecualizacion - Parte Real'),xlim([0,50]);
subplot(212),stem(imag(uk)),title('Observaciones sin ecualizacion - Parte Imaginaria'),xlim([0,50]);
thetaZF=zeros(1,length(uk));
thetaZF(1)=(uk(1+S)/valoresR(1));
thetaZF(2)=(uk(2+S)-valoresR(2)*thetaZF(1)/valoresR(1));
thetaZF(3)=(uk(3+S)-valoresR(2)*thetaZF(2)-valoresR(3)*thetaZF(1))/valoresR(1);
for i=4:length(uk)-S
    thetaZF(i)=(uk(i+S)-valoresR(2)*thetaZF(i-1)-valoresR(3)*thetaZF(i-2)/valoresR(1));
end

figure()
subplot(221),stem(real(simbolos_generados)),title('Simbolos generados - Parte Real'),xlim([0,50]);
subplot(222),stem(imag(simbolos_generados)),title('Simbolos generados - Parte Imaginaria'),xlim([0,50]);
subplot(223),stem(real(thetaZF)),title('Variable de decision ZF  - Part Real'),xlim([0,50]);
subplot(224),stem(imag(thetaZF)),title('Variable de decision ZF - Parte Imaginaria'),xlim([0,50]);

%% Decision sin ecualizacion - Distancia minima

% simbolos_recibidos=uk(S+1:S+length(simbolos_generados));
% % Inicializar variables para almacenar la distancia mínima y los índices correspondientes
% distanciaMinima = Inf;
% indicesMinimos = [];
% for j = 1:length(simbolos_recibidos)
%     for i = 1:length(constelacion)
%         distancia = abs(simbolos_recibidos(j) - constelacion(i));
% 
%         % Verificar si la distancia actual es menor que la distancia mínima
%             if distancia < distanciaMinima
%                     distanciaMinima = distancia;
%                     indicesMinimos(j) = i;
%             end
%     end
%             distanciaMinima = Inf;
% end
% 
% 
% % Mostrar los resultados
% vector_distancias=[];
% for j = 1:length(simbolos_recibidos)
%     vector_distancias(j)=constelacion(indicesMinimos(j));
% end
% %Asigna un simbolo de la constelacion a los simbolos recibidos (decision)
% simbolos_decision=vector_distancias;
% 
% %Realiza el demapeo de los simbolos recibidos en Bits
% demap= qamdemod(simbolos_decision,M,'gray');
% bits_bin=dec2bin(demap);
% bits_matrix=reshape(bits_bin.', 1, []);
% bits_demap_transpose=str2num(bits_matrix.');

%% Decision con ecualizacion ZF - Distancia minima
simbolos_recibidosZF=thetaZF(1:length(simbolos_generados));
% Inicializar variables para almacenar la distancia mínima y los índices correspondientes
distanciaMinima = Inf;
indicesMinimosZF = [];
for j = 1:length(simbolos_recibidosZF)
    for i = 1:length(constelacion)
        distancia = abs(simbolos_recibidosZF(j) - constelacion(i));

        % Verificar si la distancia actual es menor que la distancia mínima
            if distancia < distanciaMinima
                    distanciaMinima = distancia;
                    indicesMinimosZF(j) = i;
            end
    end
            distanciaMinima = Inf;
end


% Mostrar los resultados
vector_distanciasZF=[];
for j = 1:length(simbolos_recibidosZF)
    vector_distanciasZF(j)=constelacion(indicesMinimosZF(j));
end
%Asigna un simbolo de la constelacion a los simbolos recibidos (decision)
simbolos_decisionZF=vector_distanciasZF;
%Realiza el demapeo de los simbolos recibidos en Bits
demap= qamdemod(simbolos_decisionZF,M,'gray');
bits_bin=dec2bin(demap);
bits_matrix=reshape(bits_bin.', 1, []);
bits_demap_transpose=str2num(bits_matrix.');

figure();
scatter(real(simbolos_recibidosZF), imag(simbolos_recibidosZF), 'filled');
title(['Diagrama de Constelación Despues del ecualizador ',num2str(M), ' -QAM']);
xlabel('Parte Real');
ylabel('Parte Imaginaria');
grid on;
[Nerr,SER]=symerr(simbolos_generados,simbolos_decisionZF);

BER = SER*log2(M);
disp('Cantidad de Errores:')
disp(Nerr)
disp('BER:')
disp(BER)

%% recuperamos la imagen
vector_rojo_rx = bits_demap_transpose(1:columnas);
vector_verde_rx = bits_demap_transpose(columnas+1:2*columnas);
vector_azul_rx = bits_demap_transpose(2*columnas+1:end);

vector_rojo_rx_reshape = reshape(vector_rojo_rx, filas0, columnas0);
vector_verde_rx_reshape = reshape(vector_verde_rx, filas0, columnas0);
vector_azul_rx_reshape = reshape(vector_azul_rx,  filas0, columnas0);


imagen_rx=cat(3,vector_rojo_rx_reshape,vector_verde_rx_reshape,vector_azul_rx_reshape);
imagen_rx_256 = imagen_rx * 256;

figure(4)
subplot(2,1,1);
imshow(imagen_rgb);
title('Imagen Original')
subplot(2,1,2);
imshow(imagen_rx_256);
title('Imagen Reconstruida')
