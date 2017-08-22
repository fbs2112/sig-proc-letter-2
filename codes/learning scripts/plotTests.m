clear;
close all;
clc;



addpath(['.' filesep 'results']);
addpath(['..' filesep 'Utils' filesep]);

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;


figProp = struct( 'size' , 24 , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);  


fileIndex = 1;


colorCell = {'b','r','y'};

for k = 1:length(fileIndex)
    load teste2.mat;
    for i = 1:size(e3,4)
        figure
        for j = 1:size(e3,3)
            x = e3{1,1,j,i};
            blindItAux = blindIt(:,1,1,j,i);
            meanBlindIt(k,j,i) = round(mean(blindItAux(blindItAux~=0)));
            aux = find(x,1);
%             if j~=size(e3,3) && i ~= size(e3,4) && k~=2
                plot(10*log10(x(aux:end)));
%                 line([meanBlindIt(k,j,i) meanBlindIt(k,j,i)], [-20 10],'Color',colorCell{j});
                 if i > 1 && meanBlindIt(k,j,i)
                    upCountTrans(k,j,i) = mean(meanCount{1,1,j,i}(aux:meanBlindIt(k,j,i)-1))*100;
                    upCountSS(k,j,i) = mean(meanCount{1,1,j,i}(meanBlindIt(k,j,i):end))*100;
                else
                    upCount(k,j) = mean(meanCount{1,1,j,i})*100;
                end
                hold on
%             end
        end
        H = legend('$\mathrm{MI} = 0.05$','$\mathrm{MI} = 0.075$','$\mathrm{MI} = 0.1$');
        set(H,'interpreter','latex')
        
        ylabel('MSE [dB]','interpreter','latex');
        xlabel('Iterations','interpreter','latex');
        line([meanBlindIt(k,1,i) meanBlindIt(k,1,i)], [-20 10],'Color',colorCell{1});
        line([meanBlindIt(k,2,i) meanBlindIt(k,2,i)], [-20 10],'Color',colorCell{2});
        line([meanBlindIt(k,3,i) meanBlindIt(k,3,i)], [-20 10],'Color',colorCell{3});
        xlim([0 5000]);
%         formatFig( gcf ,['.' filesep 'figs' filesep '2017-08-09' filesep 'mse' num2str(i) '_' num2str(k)],'en' , figProp );

    end
end






for k = 1:length(fileIndex)
    load teste2.mat;
    for i = 1:size(e3,4)
        figure
        for j = 1:size(e3,3)
            x = e3{1,1,j,i};
            blindItAux = blindIt(:,1,1,j,i);
            meanBlindIt(k,j,i) = round(mean(blindItAux(blindItAux~=0)));
            aux = find(x,1);
%             if j~=size(e3,3) && i ~= size(e3,4) && k~=2
                plot(10*log10(x(aux:end)));
%                 line([meanBlindIt(k,j,i) meanBlindIt(k,j,i)], [-20 10],'Color',colorCell{j});
                 if i > 1 && meanBlindIt(k,j,i)
                    upCountTrans(k,j,i) = mean(meanCount{1,1,j,i}(aux:meanBlindIt(k,j,i)-1))*100;
                    upCountSS(k,j,i) = mean(meanCount{1,1,j,i}(meanBlindIt(k,j,i):end))*100;
                else
                    upCount(k,j) = mean(meanCount{1,1,j,i})*100;
                end
                hold on
%             end
        end
        H = legend('$\mathrm{MI} = 0.05$','$\mathrm{MI} = 0.075$','$\mathrm{MI} = 0.1$');
        set(H,'interpreter','latex')
        
        ylabel('MSE [dB]','interpreter','latex');
        xlabel('Iterations','interpreter','latex');
        line([meanBlindIt(k,1,i) meanBlindIt(k,1,i)], [-20 10],'Color',colorCell{1});
        line([meanBlindIt(k,2,i) meanBlindIt(k,2,i)], [-20 10],'Color',colorCell{2});
        line([meanBlindIt(k,3,i) meanBlindIt(k,3,i)], [-20 10],'Color',colorCell{3});
        xlim([0 5000]);
%         formatFig( gcf ,['.' filesep 'figs' filesep '2017-08-09' filesep 'mse' num2str(i) '_' num2str(k)],'en' , figProp );

    end
end

clear;
close all;

fileIndex = 1;


colorCell = {'b','r','y'};


for k = 1:length(fileIndex)
    load resultsDFEFF.mat;
    for i = 1:size(e3,5)
        figure
        for j = 1:size(e3,4)
            x = e3{1,1,1,j,i};
            blindItAux = blindIt(:,1,1,1,j,i);
            meanBlindIt(k,j,i) = round(mean(blindItAux(blindItAux~=0)));
            aux = find(x,1);
            %             if j~=size(e3,3) && i ~= size(e3,4) && k~=2
            plot(10*log10(x(aux:end)));
            if isnan(meanBlindIt(k,j,i))
                meanBlindIt(k,j,i) = 0;
            end
            %                 line([meanBlindIt(k,j,i) meanBlindIt(k,j,i)], [-20 10],'Color',colorCell{j});
            if i > 1 && meanBlindIt(k,j,i)
                upCountTrans(k,j,i) = mean(meanCount{1,1,1,j,i}(aux:meanBlindIt(k,j,i)-1))*100;
                upCountSS(k,j,i) = mean(meanCount{1,1,1,j,i}(meanBlindIt(k,j,i):end))*100;
            else
                upCount(k,j) = mean(meanCount{1,1,1,j,i})*100;
            end
            hold on
            %             end
        end
        H = legend('$\mathrm{MI} = 0.05$','$\mathrm{MI} = 0.075$','$\mathrm{MI} = 0.1$');
        set(H,'interpreter','latex')
        
        ylabel('MSE [dB]','interpreter','latex');
        xlabel('Iterations','interpreter','latex');
        line([meanBlindIt(k,1,i) meanBlindIt(k,1,i)], [-20 10],'Color',colorCell{1});
        line([meanBlindIt(k,2,i) meanBlindIt(k,2,i)], [-20 10],'Color',colorCell{2});
        line([meanBlindIt(k,3,i) meanBlindIt(k,3,i)], [-20 10],'Color',colorCell{3});
        xlim([0 5000]);
        %         formatFig( gcf ,['.' filesep 'figs' filesep '2017-08-09' filesep 'mse' num2str(i) '_' num2str(k)],'en' , figProp );
        
    end
end





rmpath(['.' filesep 'results']);
rmpath(['..' filesep 'Utils' filesep]);
