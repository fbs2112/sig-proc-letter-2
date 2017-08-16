clear;
close all;
clc;



addpath(['.' filesep 'results']);
addpath(['..' filesep 'Utils' filesep]);

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;


figProp = struct( 'size' , 24 , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);  

% load NLMS.mat;
% 
% it  = 4000;
% 
% 
% for i = 1:size(e5,2)
%     figure
%     
%     aux = find(e5{1,i}{1}{1},1);
%     plot(10*log10(e5{1,i}{1}{1}(aux:end)));
%     
%     c1(i) = mean(meanCount3{1,i}{1}{1}(aux:it));
% end
% 
% 
% load PNLMS.mat;
% 
% 
% for i = 1:size(e5,2)
%     figure
%     
%     aux = find(e5{1,i}{1}{1},1);
%     plot(10*log10(e5{1,i}{1}{1}(aux:end)));
%     c2(i) = mean(meanCount3{1,i}{1}{1}(aux:it));
% end
% 
% 
% 
% load DFE_PNLMS.mat;
% 
% 
% for i = 1:size(e5,1)
%     figure
%     
%     aux = find(e5{i,1}{1}{1},1);
%     plot(10*log10(e5{i,1}{1}{1}(aux:end)));
%     c3(i) = mean(meanCount3{i,1}{1}{1}(aux:it));
% end
% 
% fileIndex = [1 5];
% 
% 
% colorCell = {'b','r','y'};
% 
% for k = 1:length(fileIndex)
%     load(['results0' num2str(fileIndex(k)) '.mat']);
%     for i = 1:size(e3,4)
%         figure
%         for j = 1:size(e3,3)
%             x = e3{1,1,j,i};
%             blindItAux = blindIt(:,1,1,j,i);
%             meanBlindIt(k,j,i) = round(mean(blindItAux));
%             aux = find(x,1);
% %             if j~=size(e3,3) && i ~= size(e3,4) && k~=2
%                 plot(10*log10(x(aux:end)));
%                 line([meanBlindIt(k,j,i) meanBlindIt(k,j,i)], [-20 10],'Color',colorCell{j});
%                  if i > 1 && meanBlindIt(k,j,i)
%                     upCountTrans(k,j,i) = mean(meanCount{1,1,j,i}(aux:meanBlindIt(k,j,i)-1))*100;
%                     upCountSS(k,j,i) = mean(meanCount{1,1,j,i}(meanBlindIt(k,j,i):end))*100;
%                 else
%                     upCount(k,j) = mean(meanCount{1,1,j,i})*100;
%                 end
%                 hold on
% %             end
%         end
%         H = legend('$\mathrm{MI} = 0.05$','$\mathrm{MI} = 0.075$','$\mathrm{MI} = 0.1$');
%         set(H,'interpreter','latex')
%         
%         ylabel('MSE [dB]','interpreter','latex');
%         xlabel('Iterations','interpreter','latex');
%         formatFig( gcf ,['.' filesep 'figs' filesep '2017-08-08' filesep 'mse' num2str(i) '_' num2str(k)],'en' , figProp );
% 
%     end
% end

% clear 
figProp = struct( 'size' , 24 , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);  

% 
% colorCell = {'b','r','y'};
% 
% fileIndex = [2 3 4 6 7 8];
% 
% for k = 1:length(fileIndex)
%     load(['results0' num2str(fileIndex(k)) '.mat']);
%     for i = 1:size(e3,5)
%         figure
%         for j = 1:size(e3,4)
%             x = e3{1,1,1,j,i};
%             blindItAux = blindIt(:,1,1,1,j,i);
%             meanBlindIt(k,j,i) = round(mean(blindItAux));
%             aux = find(x,1);
%             plot(10*log10(x(aux:end)));
%             line([meanBlindIt(k,j,i) meanBlindIt(k,j,i)], [-15 10],'Color',colorCell{j});
%             if i > 1 && meanBlindIt(k,j,i)
%                 upCountTrans(j,i) = mean(meanCount{1,1,1,j,i}(aux:meanBlindIt(k,j,i)-1))*100;
%                 upCountSS(j,i) = mean(meanCount{1,1,1,j,i}(meanBlindIt(k,j,i):end))*100;
%             else
%                 upCount(j) = mean(meanCount{1,1,1,j,i})*100;
%             end
%             hold on
% 
%         end
%         
%         H = legend('$\mathrm{MI} = 0.05$','$\mathrm{MI} = 0.075$','$\mathrm{MI} = 0.1$');
%         set(H,'interpreter','latex')
%         
%         ylabel('MSE [dB]','interpreter','latex');
%         xlabel('Iterations','interpreter','latex');
% %         formatFig( gcf ,['.' filesep 'figs' filesep '2017-08-08' filesep 'mseDFE' num2str(i) '_' num2str(k)],'en' , figProp );
%         
%     end
%     
% end
% 
% 

 colorCell = {'b','r','y'};

for k = 1:1
    load teste.mat;
    for i = 1:size(e3,4)
        figure
        for j = 1:size(e3,3)
            x = e3{1,1,j,i};
            blindItAux = blindIt(:,1,1,j,i);
            meanBlindIt(j,i) = round(mean(blindItAux));
            aux = find(x,1);
%             if j~=size(e3,3) && i ~= size(e3,4) && k~=2
                plot(10*log10(x(aux:end)));
                line([meanBlindIt(j,i) meanBlindIt(j,i)], [-20 10],'Color',colorCell{j});
                 if i > 1 && meanBlindIt(j,i)
                    upCountTrans(k,j,i) = mean(meanCount{1,1,j,i}(aux:meanBlindIt(j,i)-1))*100;
                    upCountSS(k,j,i) = mean(meanCount{1,1,j,i}(meanBlindIt(j,i):end))*100;
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

    end
end








rmpath(['.' filesep 'results']);
rmpath(['..' filesep 'Utils' filesep]);
