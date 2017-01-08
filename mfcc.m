% mfcc - Mel frequency cepstrum coefficient analysis.
% % function [ceps] = mfcc(trama,samplingRate)
% Find the cepstral coefficients (ceps) corresponding to the
% input (trama) that contains the samples of a frame from a signal rated to
% samplingRate  (example: 11025)


function [ceps] = mfcc(trama,samplingRate)


%	Filter bank parameters
lowestFrequency = 133.3333;
linearFilters = 13;
linearSpacing = 66.66666666;
logFilters = 27;
logSpacing = 1.0711703;
fftSize = 512;
cepstralCoefficients = 13;
windowSize = 400;
windowSize = 256;		% Standard says 400, but 256 makes more sense
				% Really should be a function of the sample
				% rate (and the lowestFrequency) and the
				% frame rate.

                
% Keep this around for later....
totalFilters = linearFilters + logFilters;

% Now figure the band edges.  Interesting frequencies are spaced
% by linearSpacing for a while, then go logarithmic.  First figure
% all the interesting frequencies.  Lower, center, and upper band
% edges are all consequtive interesting frequencies. 

freqs = lowestFrequency + (0:linearFilters-1)*linearSpacing;
freqs(linearFilters+1:totalFilters+2) = ...
		      freqs(linearFilters) * logSpacing.^(1:logFilters+2);

lower = freqs(1:totalFilters);
center = freqs(2:totalFilters+1);
upper = freqs(3:totalFilters+2);

% We now want to combine FFT bins so that each filter has unit
% weight, assuming a triangular weighting function.  First figure
% out the height of the triangle, then we can figure out each 
% frequencies contribution
mfccFilterWeights = zeros(totalFilters,fftSize);
triangleHeight = 2./(upper-lower);
fftFreqs = (0:fftSize-1)/fftSize*samplingRate;

for chan=1:totalFilters
	mfccFilterWeights(chan,:) = (fftFreqs > lower(chan) & fftFreqs <= center(chan)).* triangleHeight(chan).*(fftFreqs-lower(chan))/(center(chan)-lower(chan)) + (fftFreqs > center(chan) & fftFreqs < upper(chan)).* triangleHeight(chan).*(upper(chan)-fftFreqs)/(upper(chan)-center(chan));
end
%semilogx(fftFreqs,mfccFilterWeights')
%axis([lower(1) upper(totalFilters) 0 max(max(mfccFilterWeights))])


% Figure out Discrete Cosine Transform.  We want a matrix
% dct(i,j) which is totalFilters x cepstralCoefficients in size.
% The i,j component is given by 
%                cos( i * (j+0.5)/totalFilters pi )
% where we have assumed that i and j start at 0.
mfccDCTMatrix = 1/sqrt(totalFilters/2)*cos((0:(cepstralCoefficients-1))' * (2*(0:(totalFilters-1))+1) * pi/2/totalFilters);
mfccDCTMatrix(1,:) = mfccDCTMatrix(1,:) * sqrt(2)/2;
%imagesc(mfccDCTMatrix);

% Filter the input with the preemphasis filter.  Also figure how
% many columns of data we will end up with.



% Allocate all the space we need for the output arrays.
ceps = zeros(cepstralCoefficients, 1);


% Invert the filter bank center frequencies.  For each FFT bin
% we want to know the exact position in the filter bank to find
% the original frequency response.  The next block of code finds the
% integer and fractional sampling positions.

% Ok, now let's do the processing.  For each chunk of data:
%    * Window the data with a hamming window,
%    * Shift it into FFT order,
%    * Find the magnitude of the fft,
%    * Convert the fft data into filter bank outputs,
%    * Find the log base 10,
%    * Find the cosine transform to reduce dimensionality.


    fftMag = abs(fft(trama,fftSize));
    t=mfccFilterWeights * fftMag;
    for i=1:size(t,1)
        for j=1:size(t,2)
            if (t(i,j) == 0)
                t(i,j) = 0.001;
            end
        end
    end

    earMag = log10(t);

    ceps = mfccDCTMatrix * earMag;

end


