clear;
clc;
close all;



addpath(['.' filesep 'LED Parameters']);

load whiteLED_334-15.mat;

%--------------Simulation parameters--------------------------------------%

maxRuns = 5000; % max runs in a single independent trial
maxIt = 1000;    %number of independent trial
signalPower = 1;    %desired input signal power
noisePower = 1e-3;  %desired measurement noise power
changingIteration = 5000;
blindIteration = 1500;
alpha = 0.9;
beta = 0.45;

%--------------Simulation parameters--------------------------------------%


%-------------------------------------------------------------------------%
%DFE set-membership

kappa = 0.5;
gamma = 1e-12;


memoryChannelLength = 3;

volterraFFFlag = 1;
volterraFBFlag = 0;

feedforwardLength = 1:5;
feedbackLength = 1:5;


feedforwardLength = 12;
feedbackLength = 12;

adaptfiltFF = zeros(length(feedforwardLength),1);
l1FF = cell(length(feedforwardLength),1);
l2FF = cell(length(feedforwardLength),1);

l1FB = cell(length(feedbackLength),1);
l2FB = cell(length(feedbackLength),1);

adaptfiltFB = zeros(length(feedbackLength),1);
adapFiltLength = zeros(length(feedforwardLength),length(feedbackLength));

for i = 1:length(feedforwardLength)
    adaptfiltFF(i) = (feedforwardLength(i)^2+feedforwardLength(i))/2 + feedforwardLength(i);
    auxMatrix = triu(ones(feedforwardLength(i)));
    [l1FF{i},l2FF{i}] = find(auxMatrix);
    for j = 1:length(feedbackLength)

        
        adaptfiltFB(j) = (feedbackLength(j)^2+feedbackLength(j))/2 + feedbackLength(j);
        

        auxMatrix = triu(ones(feedbackLength(j)));
        [l1FB{j},l2FB{j}] = find(auxMatrix);

        

        if ~volterraFFFlag
            adaptfiltFF(i) = feedforwardLength(i);
        end

        if ~volterraFBFlag
            adaptfiltFB(j) = feedbackLength(j);
        end


        adapFiltLength(i,j) = adaptfiltFF(i) + adaptfiltFB(j);
    end
end

auxMatrix = triu(ones(memoryChannelLength));
[l1Pilot,l2Pilot] = find(auxMatrix);


barGamma = 4*sqrt(5*noisePower); %threshold for set-membership purposes


numberOfBits = 2;

pamOrder = 2^numberOfBits;

SNR = db2pow(30);


%-------------------------Adaptive Filtering Parameters--------------------

delayinSamples = 13;

%-------------------------Adaptive Filtering Parameters--------------------


%-------------------------LED Parameters-----------------------------------

Poptical = @(ledLuminousEfficacy,electricalPower,k) (ledLuminousEfficacy.*electricalPower)./((1 + (ledLuminousEfficacy.*electricalPower./(maxLuminousIntensityLED/1000)).^(2*k)).^(1/(2*k)));

%-------------------------LED Parameters-----------------------------------


%-------------------------Photodiode Parameters----------------------------

A = 1e-4; %photodiode area (cm)
d = 10e-2; %distance between LED and photodiode (cm)
R = 0.5;
FOV = deg2rad(25);

%-------------------------Photodiode Parameters----------------------------


%-------------------------Pre Amplifier Parameters-------------------------

%-------------------------Pre Amplifier Parameters-------------------------


%-------------------------Transmission Parameters--------------------------

kNonLinearity = 2;

LEDfreqRespPoints = 1000;

fs = 2e6;

theta = 0;
phi = 0;

H_0 = A/d^2 * (n+1)/(2*pi) * cos(phi)^n * cos(theta) * rectangularPulse(-1,1,theta/FOV);

VDC = 3.25;
maxAbsoluteValueModulation = 3;

maxModulationIndex = (maxLEDVoltage - VDC)/VDC;
modulationIndexVector = [0.05 0.075 0.1];


save(['.' filesep 'simParameters' filesep 'paramDFE_FF.mat']);

rmpath(['.' filesep 'LED Parameters']);

