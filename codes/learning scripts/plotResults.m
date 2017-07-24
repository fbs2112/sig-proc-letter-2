clear;
close all;
clc;



addpath(['.' filesep 'results']);


load NLMS.mat;

it  = 4000;


for i = 1:size(e5,2)
    figure
    
    aux = find(e5{1,i}{1}{1},1);
    plot(10*log10(e5{1,i}{1}{1}(aux:end)));
    
    c1(i) = mean(meanCount3{1,i}{1}{1}(aux:it));
end


load PNLMS.mat;


for i = 1:size(e5,2)
    figure
    
    aux = find(e5{1,i}{1}{1},1);
    plot(10*log10(e5{1,i}{1}{1}(aux:end)));
    c2(i) = mean(meanCount3{1,i}{1}{1}(aux:it));
end



load DFE_PNLMS.mat;


for i = 1:size(e5,1)
    figure
    
    aux = find(e5{i,1}{1}{1},1);
    plot(10*log10(e5{i,1}{1}{1}(aux:end)));
    c3(i) = mean(meanCount3{i,1}{1}{1}(aux:it));
end






rmpath(['.' filesep 'results']);
