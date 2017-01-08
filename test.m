function [] = test( training_labels, training_objects, type, showKNN )
% training_labels es un vector con las labels de cada training_object.
% training_object es una matriz donde cada columna tiene una label y
% representa la feature de un trozo de audio.
%Type = 1 -> GTCC
%Type = 0 -> MFCC
% showKNN vale 0 o 1 e indica si se quiere mostrar el dibujo de KNN.
%Ya hemos entrenado, ahora pasamos a clasificar.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%.............................SETUP..............................%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%removeNoise serveix per indicar que volem quedar-nos només amb la part on
%l'ocell està cantant. Si val 1 s'eliminen els moments en que no canten els ocells. Si val 0 no s'eliminara.
removeNoise = 1;

frag = 200; %Mida dels fragments en el que dividirem l'audio en ms.
step = frag/2; %el overlap que queremos en ms. 50% en este caso.

directorioAudios = '/Users/Esterwen/Documents/MATLAB/TFG/ester/test';
list = dir([directorioAudios filesep '*.mp3']);
fileLabels = fopen('/Users/Esterwen/Documents/MATLAB/TFG/ester/labels_test.txt','r');
C = textscan(fileLabels,'%s %s');

fclose(fileLabels);
total_labels = sort(unique(training_labels));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%.............................INIT VARS..........................%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fsaux = 0; %inicialitzem la variable.
length_total_labels = length(total_labels);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%.............................TEST...............................%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(list)        
    totalFiles = length(list);
    h = waitbar(0,'Initializing waitbar...');
    for fileId=1:totalFiles
        waitbar(fileId/totalFiles,h,sprintf('%d / %d audios analizados...',fileId, totalFiles))
        [filePath,fileName,ext] = fileparts(list(fileId).name);
        for j=1:totalFiles
            
            if strcmp(fileName,C{2}{j})
                label=C{1}{j};
                
            end
        end

           if strcmp(ext,'.mp3')
                %Leemos los ficheros de audio
                [directorioAudios filesep list(fileId).name];
                [audio,fs] = audioread([directorioAudios filesep list(fileId).name]);
                
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
                %Posem a zero la variable que conta quantes labels hi ha
                %hagut de cada audio per a que es reinicii amb cada audio.
                total_labels_num = zeros (length_total_labels);

                feature_test = zeros(13,Nit);
                for j = 0: Nit-1
                    init=floor(j*stepmostres + 1);
                    fin= floor(init+ fragmostres -1);
                    %ventaneamos el fragmento de señal correspondiente.
                    audiofrag = audio(init:fin,1) .* window;
                    %calculamos mfcc o gtcc ya que los bancos ya son
                    %diferentes
                    lvl(j+1) = mean(audio(init:fin).^2);
                    if removeNoise == 1
                        if lvl(j+1)>audiolvl
                            %es un pajaro, entrenamos!
                            feature_test(:,j+1) = applyFilterBanks(audiofrag,fs,FilterWeights,DCTMatrix);   
                            [assigned_label, right_or_wrong] = KNN(training_objects , training_labels', 3, 'euclidean', feature_test(:,j+1), label, 0);
                            %contem quants n'hi ha de cada label. 
                            for k = 1 : length_total_labels
                                if strcmp(total_labels{k},assigned_label)
                                    total_labels_num(k) = total_labels_num(k) + 1;
                                end
                            end
                        end
                    else
                        %entrenamos siempre
                        feature_test(:,j+1) = applyFilterBanks(audiofrag,fs,FilterWeights,DCTMatrix);   
                        [assigned_label, right_or_wrong] = KNN(training_objects , training_labels', 3, 'euclidean', feature_test(:,j+1), label, 0);
                        %contem quants n'hi ha de cada label. 
                        for k = 1 : length_total_labels
                            if strcmp(total_labels{k},assigned_label)
                               total_labels_num(k) = total_labels_num(k) + 1;
                            end
                        end
                    end
                   
                    
                end

                
                
                
                
                fprintf('audio %d: \n', fileId);
                disp(label);
                for l = 1 : length_total_labels
                    disp([total_labels(l)' '=' total_labels_num(l)]);     
                end
                
                if showKNN ==1
                    figure,  
                    scatter(feature_test(1,:), feature_test(2,:),'+','m');
                    hold on;
                    gscatter(training_objects(1,:), training_objects(2,:),training_labels','rgbc', 'xoxo');
                end
               
           end
           
    end
    close(h);
end


end

