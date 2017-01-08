function [mfccFilterWeights,mfccDCTMatrix]= mfccFilterBanks(fs,fragmostres)

fftSize = fragmostres;


% Numero de coeficients que vols obtenir al final.
cepstralCoefficients =13;
% freq de mostratge
samplingRate = fs;


lowestFrequency=20;
linearFilters=14;
linearSpacing=200/3;
fc = lowestFrequency + (linearFilters - 1)*linearSpacing;%fins on arriben els filtres lineals
HighFreq = 20000;%El seu valor serï¿½ de 20K per samplingRate=44100. Si variem aquest valor LogFilters variarï¿½ tambï¿½.
logFilters=27;
logSpacing=exp(log(HighFreq/fc)/(logFilters+1));
CoefsN=100;
cepstral0=1;   %0 to not consider 1st cepstral (cancel DC), 1 to consider it
LogSpacing=33;

totalFilters = linearFilters+logFilters;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    FILTERBANK & DCT CONSTRUCTION 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Band edges.
freqs = lowestFrequency + (0:linearFilters+1)*linearSpacing;
% freqs(linearFilters+1:totalFilters+2) = ...
% 		      freqs(linearFilters) * logSpacing.^(1:logFilters+2);
freqs(linearFilters+3:totalFilters+2) = ...
		      freqs(linearFilters+2) * logSpacing.^(1:(logFilters));
          
lower = freqs(1:totalFilters);
center = freqs(2:totalFilters+1);
upper = freqs(3:totalFilters+2);

% Filters triangular weighting
mfccFilterWeights = zeros(totalFilters,fftSize);
triangleHeight = 2./(upper-lower);
% Si volguessim guany 1 posariem: triangleHeight = upper-(upper-1);
fftFreqs = (0:fftSize-1)/fftSize*samplingRate;

for c=1:totalFilters
	mfccFilterWeights(c,:) = ...
  (fftFreqs > lower(c) & fftFreqs <= center(c)).* ...
   triangleHeight(c).*(fftFreqs-lower(c))/(center(c)-lower(c)) + ...
  (fftFreqs > center(c) & fftFreqs < upper(c)).* ...
   triangleHeight(c).*(upper(c)-fftFreqs)/(upper(c)-center(c));
end

mfccFilterWeights=mfccFilterWeights/max(max(mfccFilterWeights));%scale max response to 0 dB
mfccFilterWeights = mfccFilterWeights';

%PLOT PER FER PROVES
%f = (0:(fftSize/2-1))*samplingRate/fftSize;

%figure (2)
%plot(f,mfccFilterWeights(1:(fftSize/2),:));axis tight;grid
%title('Resposta freqüencial filtres MFCC','Color','b')

% Discrete Cosine Transform.  
mfccDCTMatrix = 1/sqrt(totalFilters/2)*cos((0:(cepstralCoefficients-1))' * ...
				(2*(0:(totalFilters-1))+1) * pi/2/totalFilters);
mfccDCTMatrix(1,:) = mfccDCTMatrix(1,:) * sqrt(2)/2;






end

