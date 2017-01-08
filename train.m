function[training_labels, training_objects] = train(type,show)
%show indica si se quiere mostrar el spectrograma de cada audio de train.
%Type = 1 -> GTCC 
%Type = 0 -> MFCC
%retorna un vector de training_labels on cada label està assignada a un
%vector de training_objects.
%reotrna una matriu de training_objects on cada training_object és una
%columna de X coeficients GTCC o MFCC (típicament 13).



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%.............................SETUP..............................%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%removeNoise serveix per indicar que volem quedar-nos només amb la part on
%l'ocell està cantant. Si val 1 s'eliminen els moments en que no canten els ocells. Si val 0 no s'eliminara.
%SI VOLEM FER ELS PLOTS HA d'ESTAR A 0!!
removeNoise = 1;
frag = 200; %en ms.
step = frag; %el overlap que queremos en ms.

directorioAudios = '/Users/Esterwen/Documents/MATLAB/TFG/ester/trainbons';
list = dir([directorioAudios filesep '*.mp3']);
fileLabels = fopen('/Users/Esterwen/Documents/MATLAB/TFG/ester/labelsbons.txt','r');
C = textscan(fileLabels,'%s %s');
fclose(fileLabels);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%.............................INIT VARS..........................%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%serà l'index del vector de training_objects, ja que cada audio tindra X
%mfcc segons la particio de l'audio (diferent en cada cas ja que els audios
%son de mida diferent).
index = 1;
fsaux = 0; %inicialitzem la variable.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%.............................TRAIN..............................%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
        
if ~isempty(list)        
    totalFiles = length(list);
    h = waitbar(0,'Initializing waitbar...');
    for fileId=1:totalFiles
        
        waitbar(fileId/totalFiles,h,sprintf('%d / %d audios analizados...',fileId, totalFiles))
        [filePath,fileName,ext] = fileparts(list(fileId).name);
        label = 0;
        for j=1:totalFiles
            if strcmp(fileName,C{2}{j})
                label=C{1}{j};
            end
        end

           if strcmp(ext,'.mp3')
                %Leemos los ficheros de audio
                [directorioAudios filesep list(fileId).name];
                [audio,fs] = audioread([directorioAudios filesep list(fileId).name]);
                
                %calculamos el nivel medio del audio
                audiolvl = mean(audio(:,1).^2);
                %troceamos el audio en fragmentos de longitud "frag" ms.
                fragmostres = (frag / 1000) * fs ; % dividim entre 1000 ja que ho volem el ms.
                stepmostres = fs*(step/1000); %fem el mateix amb el step de la finestra per tenir el overlap.
                lengthaudio = length(audio(:,1));
                %Creamos ventana hamming.
                window = hamming(fragmostres);
                
                %Creem el banc de filtres, només ho fem si canvia la fs
                %respecte l'audio anterior.
                if fs ~= fsaux
                    fsaux = fs;
                    if type == 1
                       [FilterWeights,DCTMatrix] = gtccFilterBanks(fs,fragmostres);
                    end
                    if type == 0
                       [FilterWeights,DCTMatrix] = mfccFilterBanks(fs,fragmostres);
                    end
                end
                

                

                
                %calculamos numero iteraciones del bucle
                Nit=floor(lengthaudio/stepmostres)-ceil(frag/step);
                lvl = zeros(Nit);
                lvlin = zeros(Nit);
                for j = 0: Nit-1
                    init=floor(j*stepmostres + 1);
                    fin= floor(init+ fragmostres -1);
                    %ventaneamos el fragmento de señal correspondiente.
                    audiofrag = audio(init:fin,1) .* window;
                    %calculamos mfcc o gtcc ya que los bancos ya son
                    %diferentes.
                    %solo entrenaremos con lo que no es ruido (miramos que
                    %el señal sea más alto que la media).
                    lvl(j+1) = mean(audio(init:fin).^2);
                    if removeNoise == 1
                        if lvl(j+1)>audiolvl
                            lvlin(j+1) = 1;
                            %es un pajaro, entrenamos!
                            training_objects(:,index) = applyFilterBanks(audiofrag,fs,FilterWeights,DCTMatrix);   
                            training_labels{index} = label;
                            index = index + 1;

                        else
                            %es ruido, no entrenamos.
                            lvlin(j+1) = 0;
                        end
                    else
                         training_objects(:,index) = applyFilterBanks(audiofrag,fs,FilterWeights,DCTMatrix);   
                         training_labels{index} = label;
                         index = index + 1;
                    end
                   
                    
                    
                   
                end
                
                %Si se quiere mostrar el espectrograma...
                if (show == 1)
                    figure
                    subplot(4,1,1);
                    %Creamos ventana hamming.
                    windowSpectrogram = hamming(441);
                    %Overlap del 50%.
                    noverlap = round(length(windowSpectrogram)/2);

                    %Calculem el espectorgrama amb la propia funció del matlab.
                    spectrogram(audio(:,1),windowSpectrogram,noverlap,[],fs,'yaxis')
                    title(label)
                    set(gca,'YScale', 'log'); % flip the Y Axis so lower frequencies are at the bottom
                    %colorbar
                    xlim([0 floor(lengthaudio/(fs))-ceil(1)]);

                    subplot(4,1,2);
                    N = length(audio(:,1));
                    t = (linspace(0, N/fs - 1, N));
                    plot (t,audio(:,1));
                    xlim([0 floor(lengthaudio/(fs))-ceil(1)]);
                
                    subplot(4,1,3);
                    plot (lvl);
                    xlim([0 floor(lengthaudio/stepmostres)-ceil(frag/step)]);
                    subplot(4,1,4);
                    plot(lvlin);
                    xlim([0 floor(lengthaudio/stepmostres)-ceil(frag/step)]);
                end
                    
            end
 
    end

end

    close(h);


end