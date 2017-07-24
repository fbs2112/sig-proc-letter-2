%Semi-blind proportionate affine projection equalizer using Volterra series and VLC
%channel

clear;
clc;
close all;

addpath(['..' filesep 'simParameters' filesep]);
addpath(['..' filesep 'Utils' filesep]);

load paramEq.mat;

delayVector = N(1)+1;

maxIt = 20;

for modulationIndexLoop = 1:length(modulationIndexVector)

    for NIndex = 1:1%length(N)
        modulationIndex = modulationIndexVector(modulationIndexLoop);
        maxVoltage = VDC*(1+modulationIndex);
        deltaV = maxVoltage - VDC;

        for delay = 1:length(delayVector)



            delayVector2 = N(NIndex) + 1;
            globalLength = maxRuns + adapFiltLength(NIndex) + max(delayVector2) - 1;


            w2 = zeros(adapFiltLength(NIndex),maxRuns,maxIt);
            for index = 1:maxIt
                index


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

                w = zeros(adapFiltLength(NIndex),globalLength) + 1e-6;
                gammaAux = zeros(globalLength,1);

                for k = (adapFiltLength(NIndex) + max(delayVector2)):globalLength
                    %
                    x(:,k) = xAux(k:-1:k-N(NIndex)+1);

                    xTDLAux = zeros(length(l1{NIndex}),1);

                    for lIndex = 1:length(l1{NIndex})
                        xTDLAux(lIndex,1) = x(l1{NIndex}(lIndex),k)*(x(l2{NIndex}(lIndex),k));
                    end

                    xAP = [x(:,k);xTDLAux];


                    y(k) = w(:,k)'*xAP;

                    if k < 1000

                        d(k) = (pilot(-delayVector(delay) + k + 1));

                    else

                        %                         d(k) = pammod(pamdemod(y(k),pamOrder,0,'gray'),pamOrder,0,'gray');

                        d(k) = pamHardThreshold(y(k));

                    end

                    e(k) = d(k) - y(k);

                    %
                    gammaAux(k+1) = alpha*gammaAux(k) + (1-alpha)*sqrt(beta*w(:,k)'*w(:,k)*noisePower);

                    barGamma = sqrt(pi)*gammaAux(k+1)/2;
                    % %
                    barGamma = 4*sqrt(5*noisePower);

                    maxError = max(abs(real(e(k))),abs(imag(conj(e(k)))));

                    if maxError > barGamma
                        mu(k) = 1 - barGamma/maxError;
                        G(:,:,k) = diag(((1 - kappa*mu(k))/adapFiltLength(NIndex)) + (kappa*mu(k)*abs(w(:,k))/norm(w(:,k),1)));
                        w(:,k+1) = w(:,k) + mu(k)*G(:,:,k)*xAP*((xAP'*G(:,:,k)*xAP+gamma*eye(1))\eye(1))*conj(e(k));
                        count(k,index) = 1;
                    else
                        mu(k) = 0;
                        count(k,index) = 0;
                        w(:,k+1) = w(:,k);
                    end
                end
                wIndex(:,:,index) = conj(w(:,1:maxRuns));
                e2(:,index) = abs(e).^2;

            end
            meanCount{delay} = mean(count,2);
            w3{delay} = mean(wIndex,3);

            e3{delay} = mean(e2,2);

        end
        w4{NIndex} = w3;
        e4{NIndex} = e3;
        meanCount2{NIndex} = meanCount;

    end
    w5{modulationIndexLoop} = w4;
    e5{modulationIndexLoop} = e4;
    meanCount3{modulationIndexLoop} = meanCount2;
    figure
    plot(10*log10(e5{modulationIndexLoop}{NIndex}{1}))
    
end


save(['.' filesep 'results' filesep 'PNLMS.mat'],'w5','e5','meanCount3');

rmpath(['..' filesep 'simParameters' filesep]);
rmpath(['..' filesep 'Utils' filesep]);


