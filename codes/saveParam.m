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
%Volterra set-membership

N = 12;

l1 = cell(length(N),1);
l2 = cell(length(N),1);
adapFiltLength = zeros(length(N),1);

for i = 1:length(N)
    auxMatrix = triu(ones(N(i)));
    [l1{i},l2{i}] = find(auxMatrix);
    adapFiltLength(i) = (N(i)^2+N(i))/2 + N(i);
end

kappa = 0.5;
gamma = 1e-12;

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


save(['.' filesep 'simParameters' filesep 'paramEq.mat']);


rmpath(['.' filesep 'LED Parameters']);



