%Semi-blind proportionate affine projection equalizer using Volterra series and VLC
%channel

clear;
clc;
close all;

addpath(['..' filesep 'simParameters' filesep]);
addpath(['..' filesep 'Utils' filesep]);

load paramDFE_FF.mat;

delayVector = feedforwardLength(1)+1;

eta = 0:0.1:0.3;

e3 = cell(length(delayVector),length(feedforwardLength),length(feedbackLength),length(modulationIndexVector),length(eta));
w3 = cell(length(delayVector),length(feedforwardLength),length(feedbackLength),length(modulationIndexVector),length(eta));
meanCount = cell(length(delayVector),length(feedforwardLength),length(feedbackLength),length(modulationIndexVector),length(eta));
blindIt = zeros(maxIt,length(delayVector),length(feedforwardLength),length(feedbackLength),length(modulationIndexVector),length(eta));

% maxIt = 20;

for etaIndex = 1:length(eta)
    
    for modulationIndexLoop = 1:length(modulationIndexVector)
      
        modulationIndex = modulationIndexVector(modulationIndexLoop);
        maxVoltage = VDC*(1+modulationIndex);
        deltaV = maxVoltage - VDC;

        for FFIndex = 1:length(feedforwardLength)
            
            for FBIndex = 1:length(feedbackLength)
            
                for delay = 1:1

                    delayVector2 = feedforwardLength(FFIndex) + 1;
                    globalLength = maxRuns + adapFiltLength(FFIndex,FBIndex) + max(delayVector2) - 1;
                    
                    wIndex = zeros(adapFiltLength(FFIndex,FBIndex),globalLength,maxIt);
                    e2 = zeros(globalLength,maxIt);
                    count = zeros(globalLength,maxIt);
                    for index = 1:maxIt
                        index

                        mu = zeros(globalLength,1);
                        d = zeros(globalLength,1);                    

                        P = zeros(adapFiltLength(FFIndex,FBIndex),adapFiltLength(FFIndex,FBIndex),globalLength);
                        P(:,:,adapFiltLength(FFIndex,FBIndex) + max(delayVector2)) = eye(adapFiltLength(FFIndex,FBIndex))*1e-6;
                        sigma = zeros(globalLength,1);
                        sigma(adapFiltLength(FFIndex,FBIndex) + delayVector(delay)) = 1;
                        delta = zeros(globalLength,1);
                        lambda = zeros(globalLength,1);
                        G = zeros(globalLength,1);

                        x = zeros(feedforwardLength(FFIndex),globalLength);
                        yHat = zeros(feedbackLength(FBIndex),1);
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

                        theta = zeros(adapFiltLength(FFIndex,FBIndex),globalLength);
                        gammaAux = zeros(globalLength,1);

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
                            
                            dfeOutput = theta(:,k)'*z;

                            if k > (adapFiltLength(FFIndex,FBIndex) + max(delayVector2)) + 100
                                medianAux(k) = median(abs(delta(k-1 - 100:k-1)));
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

                            delta(k) = d(k) - dfeOutput;

                            %
                            gammaAux(k+1) = alpha*gammaAux(k) + (1-alpha)*sqrt(beta*theta(:,k)'*theta(:,k)*noisePower);

                            barGamma = sqrt(pi)*gammaAux(k+1)/2;
                            % %
                            barGamma = 4*sqrt(5*noisePower);

                            maxError = max(abs(real(delta(k))),abs(imag(conj(delta(k)))));

                            if abs(delta(k)) > barGamma
                                G(k) = z.'*P(:,:,k)*conj(z);
                                lambda(k) = (1/G(k))*((abs(delta(k))/barGamma) - 1);
                                lambda2 = 1/lambda(k);


                                P(:,:,k+1) = lambda(k)*(P(:,:,k) - (P(:,:,k)*conj(z)*z.'*P(:,:,k))/(lambda2+G(k)));


                                theta(:,k+1) = theta(:,k) + P(:,:,k+1)*conj(z)*delta(k);

                                sigma(k+1) = sigma(k) - (lambda(k)*delta(k)^2)/(1+lambda(k)*G(k)) + lambda(k)*delta(k)^2;

                                count(k,index) = 1;
                            else
                                lambda(k) = 0;
                                P(:,:,k+1) = P(:,:,k);
                                theta(:,k+1) = theta(:,k);
                                sigma(k+1) = sigma(k);
                            end
                        end
                        wIndex(:,:,index) = conj(theta(:,1:globalLength));
                        e2(:,index) = abs(delta).^2;

                    end
                    meanCount{delay,FFIndex,FBIndex,modulationIndexLoop,etaIndex} = mean(count,2);
                    w3{delay,FFIndex,FBIndex,modulationIndexLoop,etaIndex} = mean(wIndex,3);
                    
                    e3{delay,FFIndex,FBIndex,modulationIndexLoop,etaIndex} = mean(e2,2);

                end
    %             w4{NIndex} = w3;
    %             e4{NIndex} = e3;
    %             meanCount2{NIndex} = meanCount;

            end
        end
%         w5{modulationIndexLoop} = w4;
%         e5{modulationIndexLoop} = e4;
%         meanCount3{modulationIndexLoop} = meanCount2;
    end
%      w6{etaIndex} = w5;
%      e6{etaIndex} = e5;
%      meanCount4{etaIndex} = meanCount3;
end

% for i = 1:length(eta)
%     figure;
%     x = e3{1,1,1,1,i};
%     aux = find(x,1);
%     plot(10*log10(x(aux:end)))
% end
% 
% legend(eta.')
save(['.' filesep 'results' filesep 'results11.mat'],'w3','e3','meanCount','blindIt');

rmpath(['..' filesep 'simParameters' filesep]);
rmpath(['..' filesep 'Utils' filesep]);


