%This scripts generates bit error rates using pre trained adaptive filters
%as equalyzers

clear;
clc;
close all;

addpath(['..' filesep 'learning scripts' filesep 'results']);
addpath(['..' filesep 'Utils']);
addpath(['..' filesep 'berParameters']);
addpath(['..' filesep 'simParameters']);


load paramEq.mat;
load param_feedforwardEq.mat;
load results01.mat;

eta = 0:0.1:0.3;

ber = zeros(length(eta),length(modulationIndexVector),size(w3,2),length(SNR));


for etaIndex = 1:length(eta)
    
    for modulationIndexLoop = 1:length(modulationIndexVector)
        
        modulationIndex = modulationIndexVector(modulationIndexLoop);
        maxVoltage = VDC*(1+modulationIndex);
        deltaV = maxVoltage - VDC;
    
        for SNRIndex = 1:length(SNR)
            for NIndex = 1:size(w3,2)
                
                equalizerFilter = squeeze(w3(1,NIndex,modulationIndexLoop,etaIndex));
                equalizerFilter = equalizerFilter{1}(:,end);
                berAux = zeros(monteCarloLoops,1);
               
                for index = 1:monteCarloLoops
                    index
                    equalizedSignal = zeros(numberOfSymbols,1);

                    binaryInputData = randi([0,1],blockLength + 100,1);
                    binaryInputData = reshape(binaryInputData,[],numberOfBits);
                    deciInputData = bi2de(binaryInputData);    
                    pilot = pammod(deciInputData,2^numberOfBits,0,'gray');
                    convLength = length(pilot) + LEDfreqRespPoints -1;
                    NFFT = 2^nextpow2(convLength);
                    
                    pilotFreq = fft(pilot,NFFT);
                    
                    f = fs/2*linspace(0,1,NFFT/2 + 1)*2*pi;
                    
                    w = [-fliplr(f(2:end-1)) f];
                    
                    LEDResp = freqRespLED(w);
                    
                    filteredVinAux = real(ifft(pilotFreq.*fftshift(LEDResp)));
                    
                    filteredVin = filteredVinAux(1:length(pilot));
                    
                    VoltageConstant = modulationIndex*maxVoltage/((1+modulationIndex)*max(filteredVin));
                    
                    filteredVin = filteredVin*VoltageConstant + VDC;
                    
                    iLEDOutput = I_V_Fun(filteredVin,VT,nLED,ISat);
                    
                    eletricalPowerOutput = filteredVin.*iLEDOutput;
                    
                    opticalPowerOutput = Poptical(ledLuminousEfficacy,eletricalPowerOutput,kNonLinearity);
                    
                    opticalPowerOutputConvolved = opticalPowerOutput*H_0;
                    
                    n = randn(length(opticalPowerOutputConvolved),1); %noise signal
                    
                    receivedCurrentSignal = opticalPowerOutputConvolved*R*A;
                    receivedCurrentSignalAC = receivedCurrentSignal - mean(receivedCurrentSignal);
                    receivedCurrentSignalPower = receivedCurrentSignalAC'*receivedCurrentSignalAC/length(receivedCurrentSignal);
                    
                    powerNoiseAux = n'*n/(length(n));
                    powerNoise = (receivedCurrentSignalPower/db2pow(SNR(SNRIndex)));
                    n = n.*sqrt(powerNoise/powerNoiseAux);
                    
                    receivedVoltageSignalAux = (receivedCurrentSignal + n);
                    receivedVoltageSignalAux = receivedVoltageSignalAux - mean(receivedVoltageSignalAux);
                    receivedVoltageSignal =  receivedVoltageSignalAux*sqrt(var(pilot)/var(receivedVoltageSignalAux));
                    
                    
                    xAux = [zeros(N(NIndex)-1,1);receivedVoltageSignal];
                    
                    for k = N(NIndex):length(pilot)
                        
                        x = xAux(k:-1:k-N(NIndex)+1,1);
                        
                        xTDLAux = zeros(length(l1{NIndex}),1);
                        
                        for lIndex = 1:length(l1{NIndex})
                            xTDLAux(lIndex,1) = x(l1{NIndex}(lIndex))*(x(l2{NIndex}(lIndex)));
                        end
                        
                        xAP = [x;xTDLAux];

                        equalizedSignal(k,1) =  equalizerFilter'*xAP;

                    end
                    [corr,lags] = xcorr(equalizedSignal,xAux(N(NIndex):end,1));
                    [~,idx] = max(abs(corr));
                    delay = abs(lags(idx));
                    decDemodSignal = pamdemod(equalizedSignal,pamOrder,0,'gray');
                    binaryOutputData = de2bi(decDemodSignal,numberOfBits);

                    berAux(index) = sum(sum(abs(binaryOutputData(delay+1:(blockLength/2) + delay,:) - binaryInputData(1:(blockLength/2),:))))./blockLength;
                end

                ber(etaIndex,modulationIndexLoop,NIndex,SNRIndex) = mean(berAux);
            end

        end
        
    end
    
end


save(['.' filesep 'results' filesep 'resultsBER01.mat'],'SNR','ber');

rmpath(['..' filesep 'berParameters']);
rmpath(['..' filesep 'Utils']);
rmpath(['..' filesep 'learning scripts' filesep 'results']);
rmpath(['..' filesep 'simParameters']);


