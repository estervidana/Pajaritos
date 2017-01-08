function [ceps]= SN_MFCC_48Banks(input,fs,ERBFilterWeights,mfccDCTMatrix)

DCTWanted=1;
cepstral0Wanted = 1;

if size(input,2)<size(input,1)
    input=input';
end


% Numero de coeficients que vols obtenir al final.
cepstralCoefficients =13;


% Allocate all the space we need for the output arrays.
features= zeros(cepstralCoefficients);


w = input;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         PROCESSING 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flag=0;
for K=0:Nit-1
    init=floor(K*step*samplingRate + 1);
    fin= floor(init+ integ*samplingRate -1);
    long=fin-init+1;
    
    fftData = zeros(1,fftSize);       
        
    fftData(1:long)=w(init:fin).*hamming(long)';

    fftMag = abs(fft(fftData));
    earMag = mfccFilterWeights * fftMag';
    earMag2 = log10(earMag);
    ceps(:,K+1) = mfccDCTMatrix * earMag2;
    earMag3(:,K+1) = earMag2;
    flag=1;            

end
if flag==0
    earMag3=0;
end

return





