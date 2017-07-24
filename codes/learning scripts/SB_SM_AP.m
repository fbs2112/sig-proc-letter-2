%Semi-blind affinee projection equalizer



clear;
clc;
close all;

maxRuns = 1000;
maxIt = 20;

N = 20;

adapFiltLength = N;

alpha = 0.9;
beta = 0.45;

signalPower = 1;

SNRAux = db2pow(30);

noisePower = signalPower/SNRAux;
gamma = 1e-12;

h = [1 0.1294-0.4830*1i];

delayVector = round((N + length(h))/2);

M = 4;

for delay = 1:length(delayVector)

    for L = 0:0
%             count = zeros(maxIt,1);

        u = zeros(L+1,1);
        u(1) = 1;


        w2 = zeros(adapFiltLength,maxRuns,maxIt);
        for j = 1:maxIt
                j


            input = randi([0,M-1],maxRuns*2,1);
            pilot = qammod(input,M,0,'gray');

            pilot2 = pilot.*sqrt(signalPower/var(pilot));

            xAux2 = filter(h,1,pilot2);

            n = randn(maxRuns*2,1) + randn(maxRuns*2,1)*1i;
            powerSignal = xAux2'*xAux2./(maxRuns*2);
            powerNoiseAux = n'*n/(maxRuns*2);
            powerNoise = (powerSignal/SNRAux);
            n = n.*sqrt(powerNoise/powerNoiseAux);

            xAux = xAux2 + n;
% 
% 
            xTDL = (buffer(pilot,N,N-1,'nodelay'));

            w = zeros(adapFiltLength,maxRuns);

%                 w = zeros(N,maxRuns);
            gammaAux = zeros(maxRuns,1);

            for k = (adapFiltLength + delayVector(delay) + L + 10):maxRuns
%                 




                xAP = zeros(N,L+1);

                for l = 0:L
                    xAP(:,l+1) = xAux(k-l:-1:k-N+1-l);
                end

%                     xAP = flipud(xAux);

%                 xTDLConc = zeros(adapFiltLength,L+1);
% 
%                 for l3 = 1:L+1
%                     xTDLAux = zeros((N*N+N)/2,1);
% 
% 
% 
%                     for lIndex = 1:length(l1)
%                         xTDLAux(lIndex,1) = xAP(l1(lIndex),l3)*(xAP(l2(lIndex),l3));
% 
% %                                 xTDLAux(counterAux,1) = xTDL(l1,k+l3-1)*xTDL(l2,k+l3-1);
%                     end
% 
% % barGammaVector = 1:10;
% 
%                     xTDLConc(:,l3) = [xAP(:,l3);xTDLAux];
%     %                 xTDLConc(:,l3+k-L-1) = [xTDL(:,k+l3-1);xTDLAux];
%                 end

                xTDLConc = xAP;

%                 xAP = xTDLConc(:,k:-1:k-L);


                
                y(k) = w(:,k)'*xTDLConc(:,1);
                
                if k < 220
                
                    d(k) = (pilot(-delayVector(delay) + k + 1)); 
                    
                else
                    
%                     d(k) = qammod(qamdemod(y((-delayVector(delay) + k + 1)),64,0,'gray'),64,0,'gray');
                    
                    d(k) = qammod(qamdemod(y(k),M,0,'gray'),M,0,'gray');
                    
                end
                
                 e(k) = d(k) - y(k);
    

                gammaAux(k+1) = alpha*gammaAux(k) + (1-alpha)*sqrt(beta*w(:,k)'*w(:,k)*noisePower);

                barGamma = sqrt(pi)*gammaAux(k+1)/2;
%                 
                barGamma = 4*sqrt(5*noisePower);
                
                maxError = max(abs(real(e(k))),abs(imag(conj(e(k)))));

                if maxError > barGamma
                    mu(k) = 1 - barGamma/maxError;
%                     mu(k) = 0.2;
                    w(:,k+1) = w(:,k) + mu(k)*xTDLConc*((xTDLConc'*xTDLConc+gamma*eye(L+1))\eye(L+1))*conj(e(k))*u;
                    count(k,j) = 1;
                else
                    mu(k) = 0;
                    count(k,j) = 0;
                    w(:,k+1) = w(:,k);
                end
            end
            w2(:,:,j) = conj(w(:,1:maxRuns));
            e2(:,j) = abs(e).^2;
        end

        w3 = mean(w2,3);
        wFinal(delay,:,L+1) = w3(:,end);
        meanCount = mean(count,2);
        e3(delay,:,L+1) = mean(e2,2);

    end



end
figure
plot(10*log10(e3))


