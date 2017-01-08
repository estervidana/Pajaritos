function [features]= gtcc_antic(input,Setup)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             SETUP 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(input,2)<size(input,1)
    input=input';
end
% per tenir-ho en segons.
integ =Setup.General.frameLength/1000;
% solapament entre frames.
step =Setup.General.frameStep/1000;
% Numero de coeficients que vols obtenir al final.
cepstralCoefficients =Setup.SignalFeature.CoefsN;
% freq de mostratge
samplingRate = Setup.General.samplingRate;
%freq. maxima que hi haurà i la minima.
lowestFrequency=Setup.SignalFeature.LowFreq;
highestFrequency=Setup.SignalFeature.HighFreq;
% quants filtres dins del banc de filtres hi haurà.
BandsN=Setup.SignalFeature.BandsN;

% nombre de mostres que tindra el vector de la FFT. nextpow: arrodoniment a
% la potencia de 2 mes propera, ja que a la fft li va millor.
% Si jo vulgués passarli un diferent hauria de fer el temps de frame que li
% passo.
fftSize = 2^nextpow2(integ*samplingRate);
% cada banc de filtres té la mateixa mida que fft.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    FILTERBANK & DCT CONSTRUCTION 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% es fa servir la dct per comprimir la informació dels 48 filtres en 13
% coeficients. 
%aquest if crec que es per si vols tenir 13 filtres pasar del 1 al 14.
if Setup.SignalFeature.DCT==0
    cepstralCoefficients = cepstralCoefficients+1;
end

%paràmetres que serveixen per definir paràmetres dels bancs de filtres.
%Trets d'articles, diuen que donen bons resultats. Serveixen per construir
%les funcions de transferència del banc de filtres.
%%comencem a construir el banc de filtres.
EarQ = 9.26449;		
minBW = 24.7;

%per fer-ho amb una fórmula temporal o amb un altre mètode.
if Setup.SignalFeature.FormulaTemp==0
    order = 1;
else
    order_t=Setup.SignalFeature.order;
    K=1;
    fase=0;
end

%frequencies centrals del banc de filtres. (on situes cada filtre del banc)
cf= -(EarQ*minBW) + exp((1:BandsN)*(-log(highestFrequency + EarQ*minBW) + ...
    log(lowestFrequency + EarQ*minBW))/BandsN) * (highestFrequency + EarQ*minBW); % Central frequencies
%paràmetres per definir BW fde filtres, on estan situades, etc.

if Setup.SignalFeature.FormulaTemp==0
    order = 1;
    ERB = ((cf/EarQ).^order + minBW^order).^(1/order);
    B=1.019*2*pi*ERB;
else
    order_t=Setup.SignalFeature.order;
    ERB = ((cf/EarQ).^order_t+ minBW^order_t).^(1/order_t);
    B=1.019*ERB;
    K=1;
    fase=0;
end

T = 1/samplingRate;

if Setup.SignalFeature.FormulaTemp==0
    A0 = T;
    A2 = 0;
    B0 = 1;
    B1 = -2*cos(2*cf*pi*T)./exp(B*T);
    B2 = exp(-2*B*T);

    A11 = -(2*T*cos(2*cf*pi*T)./exp(B*T) + 2*sqrt(3+2^1.5)*T*sin(2*cf*pi*T)./exp(B*T))/2;
    A12 = -(2*T*cos(2*cf*pi*T)./exp(B*T) - 2*sqrt(3+2^1.5)*T*sin(2*cf*pi*T)./exp(B*T))/2;
    A13 = -(2*T*cos(2*cf*pi*T)./exp(B*T) + 2*sqrt(3-2^1.5)*T*sin(2*cf*pi*T)./exp(B*T))/2;
    A14 = -(2*T*cos(2*cf*pi*T)./exp(B*T) - 2*sqrt(3-2^1.5)*T*sin(2*cf*pi*T)./exp(B*T))/2;

%el guany de cada filtre.
    gain = abs((-2*exp(4*i*cf*pi*T)*T + 2*exp(-(B*T) + 2*i*cf*pi*T).*T.*(cos(2*cf*pi*T) - sqrt(3 - 2^(3/2))* ...
        sin(2*cf*pi*T))) .*(-2*exp(4*i*cf*pi*T)*T + 2*exp(-(B*T) + 2*i*cf*pi*T).*T.* ...
        (cos(2*cf*pi*T) + sqrt(3 - 2^(3/2)) * sin(2*cf*pi*T))).* (-2*exp(4*i*cf*pi*T)*T + ...
         2*exp(-(B*T) + 2*i*cf*pi*T).*T.* (cos(2*cf*pi*T) - sqrt(3 + 2^(3/2))*sin(2*cf*pi*T))) .* ...
       (-2*exp(4*i*cf*pi*T)*T + 2*exp(-(B*T) + 2*i*cf*pi*T).*T.* (cos(2*cf*pi*T) + sqrt(3 + 2^(3/2))*sin(2*cf*pi*T))) ./ ...
      (-2 ./ exp(2*B*T) - 2*exp(4*i*cf*pi*T) +  2*(1 + exp(4*i*cf*pi*T))./exp(B*T)).^4);       


    allfilts = ones(length(cf),1);
    A0=A0*allfilts;
    A2=A2*allfilts;
    B0=B0*allfilts;

    %Find out filter weightings
    impulse=[1 zeros(1,fftSize-1)];
    impResponse=zeros(length(impulse),BandsN);

    %per cada iteracio del bucle calcula la resposta impulsional de cada
    %filtre.
    for N=1:BandsN
        y1=filter([A0(N)/gain(N) A11(N)/gain(N) A2(N)/gain(N)], [B0(N) B1(N) B2(N)], impulse);
        y2=filter([A0(N) A12(N) A2(N)], [B0(N) B1(N) B2(N)], y1);  
        y3=filter([A0(N) A13(N) A2(N)], [B0(N) B1(N) B2(N)], y2);
        y4=filter([A0(N) A14(N) A2(N)], [B0(N) B1(N) B2(N)], y3);
        impResponse(:,BandsN+1-N) = y4;    %turn matrix, so lower indices realted to lower freqs

    end
else
  for j=1:BandsN   
   Anorm = zeros(1,fftSize);
%    for n=0:(fftSize-1)
%        t(n+1) = n*T;
%        g(n+1,j)=K*t(n+1)^(order_t-1)*exp(-2*pi*B(j)*t(n+1))*cos(2*pi*cf(j)*t(n+1)+fase);%funció resposta impulsional Gammatone
%        Anorm(n+1) = K*t(n+1)^(order_t-1)*exp(-2*pi*B(j)*t(n+1));
%    end
n = (0:(fftSize-1));
t = n*T;
g(:,j)=K*(t(n+1).^(order_t-1)).*exp(-2*pi*B(j)*t).*cos(2*pi*cf(j)*t+fase);%funció resposta impulsional Gammatone
Anorm = K*(t.^(order_t-1)).*exp(-2*pi*B(j)*t);

   if Setup.SignalFeature.TipusGuany
     g(:,j) = 2*g(:,j)/sum(Anorm);%normalització peak
   else
     g(:,j) = g(:,j)/sqrt(sum(g(:,j).^2)/2);%normalització per energia
   end
  
  end

end
%% Fins aquí es construeix el banc de filtres.
if Setup.SignalFeature.FormulaTemp==1
    impResponse=g;
end

%COMPROVACIONS GRÀFIQUES ENTRE FILTRES D'ORDRE 1 EN CASCADA O BANC DE
%FILTRES CALCULATS AMB FÒRMULA TEMPORAL:
% ERBFilterTemps=abs(fft(g));
% ERBFilterTemps(fftSize/2:fftSize,:)=0;
% ERBFilterWeights=abs(fft(impResponse));
% ERBFilterWeights(fftSize/2:fftSize,:)=0;    %remove the repeated spectra
% 
% f = (0:(fftSize/2-1))*samplingRate/fftSize;
% 
% figure(1)
% subplot(2,1,1),plot(f,ERBFilterTemps(1:(fftSize/2),:));axis tight;grid
% title('Resposta freqüencial filtres a partir de la fórmula temporal','Color','b')
% 
% subplot(2,1,2),plot(f,ERBFilterWeights(1:(fftSize/2),:));axis tight;grid
% title('Resposta freqüencial filtres métode Xavi Valero','Color','b')

%valor absolut de la resposta impulsional per tenir la resposta en
%frequencia de cada filtre. Es una matriu que cada columna son els valors
%de la resposta impulsional
%de cada filtre del banc.
ERBFilterWeights=abs(fft(impResponse));
%ens quedem només amb una meitat.
ERBFilterWeights(fftSize/2:fftSize,:)=0;
 
% Discrete Cosine Transform. podria ser gtccDCTMatrix. Es la matriu que permet comprimir la informació.
% En temps real sempre valdria el mateix.
mfccDCTMatrix = 1/sqrt(BandsN/2)*cos((0:(cepstralCoefficients-1))' * ...
				(2*(0:(BandsN-1))+1) * pi/2/BandsN);
mfccDCTMatrix(1,:) = mfccDCTMatrix(1,:) * sqrt(2)/2;


% figure how many columns of data we will end up with.
% Nit = numero d'iteracions que farà el filtre.
Nit=floor(length(input)/(samplingRate*step))-ceil(integ/step);

% Allocate all the space we need for the output arrays.
ceps= zeros(cepstralCoefficients, Nit);
earMag2= zeros(BandsN, Nit);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         PROCESSING 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%no cal fer for si només vols fer-ho un cop.
for K=0:Nit-1
%mostra inicial i final del fitxer input.
    init=floor(K*step*samplingRate + 1);
    fin= floor(init+ integ*samplingRate -1);
    %ester init = 1;
    %ester fin = length(input);
    long=fin-init+1;
    
    %la primera linea que hauria de pillar.
    fftData = zeros(1,fftSize);      
    %Multipliquem per la hamming per a tenir-ho enfinestrat, encara no és
    %la FFT.
    fftData(1:long)=input(init:fin).*hamming(long)';
%     N=long;
%     hamming=0.53836-0.46164*cos(2*pi*(0:(N-1))/(N-1));

   %Aqui ja fem la FFT, ens quedem amb la part real.
    fftMag=abs(fft(fftData));
    %lo que s'assembla a la oida humana. log10 per passar-ho a dB i
    %que s'assembli a la oida humana. (la multiplicacio és filtrar-ho).
    earMag = log10(ERBFilterWeights'*fftMag');
    %earMag=earMag./sqrt(earMag'*earMag); %spectra normalization
    earMag2(:,K+1)=earMag; %Test. Basura, sembla q no serveix.
    %multipliquem la matriu amb els vectors de gtcc per els coeficients de la DCT
    %per tenir el nombre de coef. que voliem.
    ceps(:,K+1) = mfccDCTMatrix * earMag;
end

if Setup.SignalFeature.DCT==0
    features=earMag2;
else
    features=ceps;
    %si es 0 vol dir que no agafa el primer coeficient, ja que el primer
    %coeficient fa més alusió al nivell de continua i els altres
    %coeficients fan alusio a tot el contingut frequencial.
    if Setup.SignalFeature.cepstral0==0
        features=features(2:(cepstralCoefficients),:);
    end

    %ESTER: Fem el promig de les features per tenir una columna amb els 13
    %gtcc
features = sum(features,2)/size(features,1);
end