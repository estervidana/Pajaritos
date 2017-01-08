function [features]= gtcc(input,fs,FilterWeights,DCTMatrix)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             SETUP 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DCTWanted=1;
cepstral0Wanted = 1;

if size(input,2)<size(input,1)
    input=input';
end


% Numero de coeficients que vols obtenir al final.
cepstralCoefficients =13;


% Allocate all the space we need for the output arrays.
features= zeros(cepstralCoefficients);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         PROCESSING 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %JA LI PASSEM EL SENYAL ENFINESTRAT, NO CAL TORNAR A APLICAR LA
    %FINESTRA.
    %Declarem la llargada de la fft.
    %fftData = zeros(1,fftSize);      
    %Multipliquem per la hamming per a tenir-ho enfinestrat, encara no és
    %la FFT.
    %fftData(1:fftSize)=input.*hamming(fftSize)';

   %Aqui ja fem la FFT, ens quedem amb la part real.
    fftMag=abs(fft(input));
    %lo que s'assembla a la oida humana. log10 per passar-ho a dB i
    %que s'assembli a la oida humana. (la multiplicacio és filtrar-ho).
    
    earMag = log10(FilterWeights'*fftMag');
    
    
    
    if DCTWanted==0
        features=earMag;
    else
    
    %multipliquem la matriu amb els vectors de gtcc per els coeficients de la DCT
    %per tenir el nombre de coef. que voliem.
    
    features = DCTMatrix * earMag;
    end

    %si es 0 vol dir que no agafa el primer coeficient, ja que el primer
    %coeficient fa més alusió al nivell de continua i els altres
    %coeficients fan alusio a tot el contingut frequencial.
    if cepstral0Wanted==0
        features=features(2:(cepstralCoefficients),:);
    end

end