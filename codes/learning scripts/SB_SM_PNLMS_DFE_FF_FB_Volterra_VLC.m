%Semi-blind proportionate affine projection equalizer using Volterra series and VLC
%channel

clear;
clc;
close all;

addpath(['..' filesep 'simParameters' filesep]);
addpath(['..' filesep 'Utils' filesep]);

load paramDFE_FF_FB.mat;

delayVector = feedforwardLength(1)+1;
eta = 0:0.1:0.3;

% maxIt = 20;

e3 = cell(length(delayVector),length(feedforwardLength),length(feedbackLength),length(modulationIndexVector),length(eta));
w3 = cell(length(delayVector),length(feedforwardLength),length(feedbackLength),length(modulationIndexVector),length(eta));
meanCount = cell(length(delayVector),length(feedforwardLength),length(feedbackLength),length(modulationIndexVector),length(eta));
blindIt = zeros(maxIt,length(delayVector),length(feedforwardLength),length(feedbackLength),length(modulationIndexVector),length(eta));

for etaIndex = 1:length(eta)
    
    for modulationIndexLoop = 1:length(modulationIndexVector)
        
        modulationIndex = modulationIndexVector(modulationIndexLoop);
        maxVoltage = VDC*(1+modulationIndex);
        deltaV = maxVoltage - VDC;
        
        for FFIndex = 1:length(feedforwardLength)
            
            for FBIndex = 1:length(feedbackLength)
                
                for delay = 1:length(delayVector)
                    
                    delayVector2 = feedforwardLength(FFIndex) + 1;
                    globalLength = maxRuns + adapFiltLength(FFIndex,FBIndex) + max(delayVector2) - 1;
                    
                    wIndex = zeros(adapFiltLength(FFIndex,FBIndex),globalLength,maxIt);
                    e2 = zeros(globalLength,maxIt);
                    count = zeros(globalLength,maxIt);
                    
                    for index = 1:maxIt
                        index
                        
                        mu = zeros(globalLength,1);
                        d = zeros(globalLength,1);
                        e = zeros(globalLength,1);
                        G = zeros(adapFiltLength(FFIndex,FBIndex),adapFiltLength(FFIndex,FBIndex),globalLength);
                        
                        x = zeros(feedforwardLength(FFIndex),globalLength);
                        yHat = zeros(feedbackLength(FBIndex),1);
                        gammaAux = zeros(globalLength,1);
                        medianAux = zeros(globalLength,1);
                        
                        input = randi([0,pamOrder-1],globalLength,1);
                        pilot = real(pammod(input,pamOrder,0,'gray'));
                        
                        %             pilot2 = pilot.*sqrt(signalPower/var(pilot));
                        
                        
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
                        powerNoise = (receivedCurrentSignalPower/SNR);
                        n = n.*sqrt(powerNoise/powerNoiseAux);
                        
                        receivedVoltageSignalAux = (receivedCurrentSignal + n);
                        receivedVoltageSignalAux = receivedVoltageSignalAux - mean(receivedVoltageSignalAux);
                        receivedVoltageSignal =  receivedVoltageSignalAux*sqrt(var(pilot)/var(receivedVoltageSignalAux));
                        
                        
                        xAux = [zeros(feedforwardLength(FFIndex)-1,1);receivedVoltageSignal];
                        
                        w = zeros(adapFiltLength(FFIndex,FBIndex) ,globalLength) + 1e-6;
                        blindFlag = 0;
                        
                        for k = (adapFiltLength(FFIndex,FBIndex) + max(delayVector2)):globalLength
                            %
                            x(:,k) = xAux(k:-1:k-feedforwardLength(FFIndex)+1);
                            
                            yHat(:,k) = (pilot(-delayVector2 + k + 1 -1:-1:-delayVector2 + k + 1 - feedbackLength(FBIndex) - 1 + 1));
                            
                            if volterraFFFlag
                                
                                aux = zeros((feedforwardLength(FFIndex)^2+feedforwardLength(FFIndex))/2,1);
                                
                                for lIndex = 1:length(l1FF{FFIndex})
                                    aux(lIndex,1) = x(l1FF{FFIndex}(lIndex),k)*(x(l2FF{FFIndex}(lIndex),k));
                                end
                                xConc = [x(:,k);aux];
                            else
                                xConc = x(:,k);
                            end
                            
                            
                            if volterraFBFlag
                                aux = zeros((feedbackLength(FBIndex)^2+feedbackLength(FBIndex))/2,1);
                                for lIndex = 1:length(l1FB{FBIndex})
                                    aux(lIndex,1) = yHat(l1FB{FBIndex}(lIndex),k)*(yHat(l2FB{FBIndex}(lIndex),k));
                                end
                                
                                yHatConc = [yHat(:,k);aux];
                            else
                                yHatConc = yHat(:,k);
                            end
                            
                            if ~volterraFFFlag && ~volterraFBFlag
                                xConc = x(:,k);
                                yHatConc = yHat(:,k);
                            end
                            
                            z = [xConc;yHatConc];
                            
                            dfeOutput = w(:,k)'*z;
                            
                            if k > (adapFiltLength(FFIndex,FBIndex) + max(delayVector2)) + 100
                                medianAux(k) = median(abs(e(k-1 - 100:k-1)));
                                if medianAux(k) <= 2*eta(etaIndex) || blindFlag == 1
                                    d(k) = pamHardThreshold(dfeOutput);
                                    
                                    if ~blindFlag
                                        blindIt(index,delay,FFIndex,FBIndex,modulationIndexLoop,etaIndex) = k;
                                    end
                                    blindFlag = 1;
                                    
                                else
                                    d(k) = (pilot(-delayVector(delay) + k + 1));
                                end
                                
                            else
                                d(k) = (pilot(-delayVector(delay) + k + 1));
                                %                         d(k) = pammod(pamdemod(y(k),pamOrder,0,'gray'),pamOrder,0,'gray');
                            end
                            
                            e(k) = d(k) - dfeOutput;
                            
                            maxError = max(abs(real(e(k))),abs(imag(conj(e(k)))));
                            
                            gammaAux(k+1) = alpha*gammaAux(k) + (1-alpha)*sqrt(beta*w(:,k)'*w(:,k)*noisePower);
                            
                            barGamma = gammaAux(k+1);
%                             barGamma = 4*sqrt(5*noisePower);
                            
                            if maxError > barGamma
                                mu(k) = 1 - barGamma/maxError;
                                G(:,:,k) = diag(((1 - kappa*mu(k))/adapFiltLength(FFIndex,FBIndex)) + (kappa*mu(k)*abs(w(:,k))/norm(w(:,k),1)));
                                w(:,k+1) = w(:,k) + mu(k)*G(:,:,k)*z*((z'*G(:,:,k)*z+gamma*eye(1))\eye(1))*conj(e(k));
                                count(k,index) = 1;
                            else
                                mu(k) = 0;
                                w(:,k+1) = w(:,k);
                                G(:,:,k) = eye(adapFiltLength(FFIndex,FBIndex));
                            end
                            
                        end
                        wIndex(:,:,index) = conj(w(:,1:globalLength));
                        e2(:,index) = abs(e).^2;
                        
                    end
                    meanCount{delay,FFIndex,FBIndex,modulationIndexLoop,etaIndex} = mean(count,2);
                    w3{delay,FFIndex,FBIndex,modulationIndexLoop,etaIndex} = mean(wIndex,3);
                    
                    e3{delay,FFIndex,FBIndex,modulationIndexLoop,etaIndex} = mean(e2,2);
                    
                end
            end
            
        end
        
    end
end

% for i = 1:length(eta)
%     plot(10*log10(e3{1,1,1,1,i}))
%     hold on
% end

% legend(eta.')


save(['.' filesep 'results' filesep 'results14.mat'],'w3','e3','meanCount','blindIt');

rmpath(['..' filesep 'simParameters' filesep]);
rmpath(['..' filesep 'Utils' filesep]);


