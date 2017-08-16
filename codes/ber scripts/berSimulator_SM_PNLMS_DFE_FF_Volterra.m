%This scripts generates bit error rates using pre trained adaptive filters
%as equalyzers

clear;
clc;
close all;

addpath(['..' filesep 'learning scripts' filesep 'results']);
addpath(['..' filesep 'Utils']);
addpath(['..' filesep 'berParameters']);
addpath(['..' filesep 'simParameters']);


load paramDFE_FF.mat;
load param_feedforwardEq.mat;
load results03.mat;

eta = 0:0.1:0.3;

ber = zeros(length(eta),length(modulationIndexVector),size(w3,2),size(w3,3),length(SNR));


for etaIndex = 1:length(eta)
    
    for modulationIndexLoop = 1:length(modulationIndexVector)
        
        modulationIndex = modulationIndexVector(modulationIndexLoop);
        maxVoltage = VDC*(1+modulationIndex);
        deltaV = maxVoltage - VDC;
        
        for SNRIndex = 1:length(SNR)
            for FFIndex = 1:length(feedforwardLength)
                for FBIndex = 1:length(feedbackLength)
                    
                    equalizerFilter = squeeze(w3(1,FFIndex,FBIndex,modulationIndexLoop,etaIndex));
                    equalizerFilter = equalizerFilter{1}(:,end);
                    berAux = zeros(monteCarloLoops,1);
                    
                    for index = 1:monteCarloLoops
                        index
                        equalizedSignal = zeros(numberOfSymbols,1);
                        
                        binaryInputData = randi([0,1],blockLength + 100,1);
                        binaryInputData = reshape(binaryInputData,[],numberOfBits);
                        deciInputData = bi2de(binaryInputData);
                        pilot = real(pammod(deciInputData,2^numberOfBits,0,'gray'));
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
                        
                        
                        xAux = [zeros(feedforwardLength(FFIndex)-1,1);receivedVoltageSignal];
                        outputFF = zeros(length(pilot),1);
                        outputFB = zeros(length(pilot),1);
                        
                        for k = feedforwardLength(FFIndex):length(pilot)
                            
                            
                            x = xAux(k:-1:k-feedforwardLength(FFIndex)+1,1);
                            
                            if volterraFFFlag
                                
                                aux = zeros((feedforwardLength(FFIndex)^2+feedforwardLength(FFIndex))/2,1);
                                
                                for lIndex = 1:length(l1FF{FFIndex})
                                    aux(lIndex,1) = x(l1FF{FFIndex}(lIndex),1)*(x(l2FF{FFIndex}(lIndex),1));
                                end
                                xConc = [x(:,1);aux];
                            else
                                xConc = x(:,1);
                            end
                            outputFF(k) =  (equalizerFilter(1:adaptfiltFF(FFIndex)))'*xConc;
                            equalizedSignal(k,1) = pamHardThreshold(outputFF(k) + outputFB(k-1));
                            
                            inputFB = equalizedSignal(k:-1:k-feedbackLength(FBIndex) + 1);
                            
                            if volterraFBFlag
                                aux = zeros((feedbackLength(FBIndex)^2+feedbackLength(FBIndex))/2,1);
                                for lIndex = 1:length(l1FB{FBIndex})
                                    aux(lIndex,1) = inputFB(l1FB{FBIndex}(lIndex),1)*(inputFB(l2FB{FBIndex}(lIndex),1));
                                end
                                
                                yHatConc = [inputFB(:,1);aux];
                            else
                                yHatConc = inputFB(:,1);
                            end
                            
                            if ~volterraFFFlag && ~volterraFBFlag
                                xConc = x(:,1);
                                yHatConc = inputFB;
                            end
                            outputFB(k) = (equalizerFilter(adaptfiltFF(FFIndex)+1:end,1))'*yHatConc;
                        end
                        
                        [corr,lags] = xcorr(equalizedSignal,xAux(feedforwardLength(FFIndex):end,1));
                        [~,idx] = max(abs(corr));
                        delay = abs(lags(idx));
                        decDemodSignal = pamdemod(equalizedSignal,pamOrder,0,'gray');
                        binaryOutputData = de2bi(decDemodSignal,numberOfBits);
                        
                        berAux(index) = sum(sum(abs(binaryOutputData(delay+1:(blockLength/2) + delay,:) - binaryInputData(1:(blockLength/2),:))))./blockLength;
                    end
                    
                    ber(etaIndex,modulationIndexLoop,FFIndex,FBIndex,SNRIndex) = mean(berAux);
                end
                
            end
        end
        
    end
    
end


save(['.' filesep 'results' filesep 'resultsBER03.mat'],'SNR','ber');

rmpath(['..' filesep 'berParameters']);
rmpath(['..' filesep 'Utils']);
rmpath(['..' filesep 'learning scripts' filesep 'results']);
rmpath(['..' filesep 'simParameters']);


