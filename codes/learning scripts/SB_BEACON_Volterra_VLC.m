%Semi-blind proportionate affine projection equalizer using Volterra series and VLC
%channel

clear;
clc;
close all;

addpath(['..' filesep 'simParameters' filesep]);
addpath(['..' filesep 'Utils' filesep]);

load paramEq.mat;

delayVector = N(1)+1;

eta = 0:0.1:0.3;

e3 = cell(length(delayVector),length(N),length(modulationIndexVector),length(eta));
w3 = cell(length(delayVector),length(N),length(modulationIndexVector),length(eta));
meanCount = cell(length(delayVector),length(N),length(modulationIndexVector),length(eta));
blindIt = zeros(maxIt,length(delayVector),length(N),length(modulationIndexVector),length(eta));

maxIt = 20;

for etaIndex = 1:1%length(eta)
    
    for modulationIndexLoop = 1:1%length(modulationIndexVector)
      
        modulationIndex = modulationIndexVector(modulationIndexLoop);
        maxVoltage = VDC*(1+modulationIndex);
        deltaV = maxVoltage - VDC;

        for NIndex = 1:1%length(N)
            
            for delay = 1:1

                delayVector2 = N(NIndex) + 1;
                globalLength = maxRuns + adapFiltLength(NIndex) + max(delayVector2) - 1;

                e2 = zeros(globalLength,maxIt);
                wIndex = zeros(adapFiltLength(NIndex),globalLength,maxIt);
                count = zeros(globalLength,maxIt);
                for index = 1:maxIt
                    index
                    
                    mu = zeros(globalLength,1);
                    d = zeros(globalLength,1);
                    e = zeros(globalLength,1);
                    
                    
                    P = zeros(adapFiltLength(NIndex),adapFiltLength(NIndex),globalLength);
                    P(:,:,adapFiltLength(NIndex) + max(delayVector2)) = eye(adapFiltLength(NIndex))*1e-6;
                    sigma = zeros(globalLength,1);
                    sigma(adapFiltLength(NIndex) + delayVector(delay)) = 1;
                    delta = zeros(globalLength,1);
                    lambda = zeros(globalLength,1);
                    G = zeros(globalLength,1);
                    
                    x = zeros(N(NIndex),globalLength);
                    medianAux = zeros(globalLength,1);
                    y = zeros(globalLength,1);
                    
                    
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


                    xAux = [zeros(N(NIndex)-1,1);receivedVoltageSignal];

                    theta = zeros(adapFiltLength(NIndex),globalLength);
                    gammaAux = zeros(globalLength,1);

                    blindFlag = 0;

                    for k = (adapFiltLength(NIndex) + max(delayVector2)):globalLength
                        %
                        x(:,k) = xAux(k:-1:k-N(NIndex)+1);

                        xTDLAux = zeros(length(l1{NIndex}),1);

                        for lIndex = 1:length(l1{NIndex})
                            xTDLAux(lIndex,1) = x(l1{NIndex}(lIndex),k)*(x(l2{NIndex}(lIndex),k));
                        end

                        xAP = [x(:,k);xTDLAux];


                        y(k) = theta(:,k)'*xAP;


                        if k > (adapFiltLength(NIndex) + max(delayVector2)) + 100
                            medianAux(k) = median(abs(e(k-1 - 100:k-1)));
                            if medianAux(k) <= 2*eta(etaIndex) || blindFlag == 1
                                d(k) = pamHardThreshold(y(k));

                                if ~blindFlag
                                    blindIt(index,delay,NIndex,modulationIndexLoop,etaIndex) = k;
                                end
                                blindFlag = 1;
                                
                            else
                                d(k) = (pilot(-delayVector(delay) + k + 1));
                            end

                        else
                            d(k) = (pilot(-delayVector(delay) + k + 1));
                            %                         d(k) = pammod(pamdemod(y(k),pamOrder,0,'gray'),pamOrder,0,'gray');
                        end

                        delta(k) = d(k) - theta(:,k)'*xAP;

                        %
                        gammaAux(k+1) = alpha*gammaAux(k) + (1-alpha)*sqrt(beta*theta(:,k)'*theta(:,k)*noisePower);

                        barGamma = sqrt(pi)*gammaAux(k+1)/2;
                        % %
                        barGamma = 4*sqrt(5*noisePower);

                        maxError = max(abs(real(delta(k))),abs(imag(conj(delta(k)))));

                        if abs(delta(k)) > barGamma
                            G(k) = xAP.'*P(:,:,k)*conj(xAP);
                            lambda(k) = (1/G(k))*((abs(delta(k))/barGamma) - 1);
                            lambda2 = 1/lambda(k);
                            
                            
                            P(:,:,k+1) = lambda(k)*(P(:,:,k) - (P(:,:,k)*conj(xAP)*xAP.'*P(:,:,k))/(lambda2+G(k)));
                            
                            
                            theta(:,k+1) = theta(:,k) + P(:,:,k+1)*conj(xAP)*delta(k);
                            
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
                meanCount{delay,NIndex,modulationIndexLoop,etaIndex} = mean(count,2);
                w3{delay,NIndex,modulationIndexLoop,etaIndex} = mean(wIndex,3);

                e3{delay,NIndex,modulationIndexLoop,etaIndex} = mean(e2,2);

            end
%             w4{NIndex} = w3;
%             e4{NIndex} = e3;
%             meanCount2{NIndex} = meanCount;

        end
%         w5{modulationIndexLoop} = w4;
%         e5{modulationIndexLoop} = e4;
%         meanCount3{modulationIndexLoop} = meanCount2;
    end
%      w6{etaIndex} = w5;
%      e6{etaIndex} = e5;
%      meanCount4{etaIndex} = meanCount3;
end

for i = 1:length(eta)
    plot(10*log10(e3{1,1,1,i}))
    hold on
end
% 
% legend(eta.')
% save(['.' filesep 'results' filesep 'results01.mat'],'w3','e3','meanCount','blindIt');

rmpath(['..' filesep 'simParameters' filesep]);
rmpath(['..' filesep 'Utils' filesep]);


