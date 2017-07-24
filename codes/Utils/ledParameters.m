%This script saves a .mat file that contains some LED parameters, which is
%used by the simulator.

clear;
clc;
close all

addpath(['.' filesep 'LED Parameters']);

load whiteLED_334-15_Param.mat;

maxLEDVoltage = 3.6; %V
minLEDVoltage = 3;
maxLEDCurrent = 0.03; %mA
minLEDCurrent = 0.004; %mA

maxElectricalPower = maxLEDVoltage*maxLEDCurrent;
minElectricalPower = minLEDCurrent*minLEDVoltage;

ISat = ISat; %saturation current
VB = 2.6; %minimum voltage for current flow 
nLED = n; %LED ideality factor
VT = 0.025; %Thermal voltage

halfAngleLED = deg2rad(15); %LED half power angle in degrees
n = -log(2)/log(cos(halfAngleLED));
luminousIntensityLED = 21375; %milicandela
maxLuminousIntensityLED = 28500;%milicandela

maxCd = 28.5; %candela
minCd = 14.25; %candela

ledLuminousEfficacy = (maxCd - minCd)/(maxElectricalPower - minElectricalPower) ; %LED luminous efficacy in cd/W


save(['.' filesep 'LED Parameters' filesep 'whiteLED_334-15.mat']);

rmpath(['.' filesep 'LED Parameters']);
