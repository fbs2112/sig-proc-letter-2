%Semi-blind affine projection equalizer using Volterra series



clear;
clc;
close all;

addpath(['..' filesep 'simParameters' filesep]);

load paramEq.mat;

delayVector = N(5)+1;

maxIt = 20;


alpha = 0.9;
beta = 0.45;

for NIndex = 5:5%length(N)


    for delay = 1:length(delayVector)
        
        

        delayVector2 = [N(NIndex)+1 N(NIndex)-2];
        globalLength = maxRuns + adapFiltLength(NIndex) + max(delayVector2) - 1;


            w2 = zeros(adapFiltLength(NIndex),maxRuns,maxIt);
            for index = 1:maxIt
                 index


                input = randi([0,pamOrder-1],globalLength,1);
                pilot = pammod(input,pamOrder,0,'gray');

                pilot2 = pilot.*sqrt(signalPower/var(pilot));



                xAux = zeros(length(pilot),size(h,2));

                for channelIndex = 1:size(h,2)
                    aux2 = zeros(length(l1Pilot),1);
                    xAux2 = zeros(length(pilot),1);

                    for i = memoryChannelLength:length(pilot) %Channel 1
                       xPilot = (pilot2(i:-1:i-memoryChannelLength+1));
                       for lIndex = 1:length(l1Pilot)
                          aux2(lIndex,1) = xPilot(l1Pilot(lIndex),1)*(xPilot(l2Pilot(lIndex),1));
                       end
                       xConc = [xPilot;(aux2)];
                       xAux2(i,1) = xConc.'*h(:,channelIndex);
                    end


        %         n = randn(globalLength,1) + randn(globalLength,1)*1i;

                    n = randn(globalLength,1);
                    powerSignal = xAux2'*xAux2./(globalLength);
                    powerNoiseAux = n'*n/(globalLength);
                    powerNoise = (powerSignal/SNR);
                    n = n.*sqrt(powerNoise/powerNoiseAux);

                    xAux(:,channelIndex) = xAux2 + n;

                end

                channelIndex = 1;

                 w = zeros(adapFiltLength(NIndex),globalLength) + 1e-6;
                gammaAux = zeros(globalLength,1);

                for k = (adapFiltLength(NIndex) + max(delayVector2)):globalLength
    %                 
                    x(:,k) = xAux(k:-1:k-N(NIndex)+1,channelIndex);

                    xTDLAux = zeros(length(l1{NIndex}),1);

                    for lIndex = 1:length(l1{NIndex})
                        xTDLAux(lIndex,1) = x(l1{NIndex}(lIndex),k)*(x(l2{NIndex}(lIndex),k));
                    end

                    xAP = [x(:,k);xTDLAux];
                   

                    y(k) = w(:,k)'*xAP;

                    if k < 2000

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

%                     barGamma = 1;

                    maxError = max(abs(real(e(k))),abs(imag(conj(e(k)))));

                    if maxError > barGamma
                        mu(k) = 1 - barGamma/maxError;
    %                     mu(k) = 0.2;
                        w(:,k+1) = w(:,k) + mu(k)*xAP*((xAP'*xAP+gamma*eye(1))\eye(1))*conj(e(k));
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

figure
plot(10*log10(e4{5}{1}))


rmpath(['..' filesep 'simParameters' filesep]);




