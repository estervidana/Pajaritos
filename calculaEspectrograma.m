function [] = calculaEspectrograma(file,label)
%Esta funci�n muestra tanto el espectrograma como la evoluci�n temporal del
%se�al correspondiente a "file".
%"file" debe incluir el path/nombrefichero.extensi�n del fichero a mostrar.
%"label ser� el titulo de la figura.
%Leemos el fichero de audio
[audio,fs] = audioread(file);

%Si l'audio es estereo el convertim en mono agafant nomes el canal 1.
audio = audio(:,1);
%Creamos ventana hamming.
window = hamming(441);
sizeWin = size(window);
%Overlap del 50%.
noverlap = round(length(window)/2);

figure
subplot(2,1,1);
%Calculem el espectorgrama amb la propia funci� del matlab.
spectrogram(audio,window,noverlap,[],fs,'yaxis')
title(label)
set(gca,'YScale', 'log'); % flip the Y Axis so lower frequencies are at the bottom
colorbar()

subplot(2,1,2);
N = length(audio);
plot ((linspace(0, N/fs, N)),audio); 

                

end

