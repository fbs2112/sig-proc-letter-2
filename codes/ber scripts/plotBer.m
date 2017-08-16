clear;
clc;
close all;


addpath(['.' filesep 'results']);

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;

figProp = struct( 'size' , 24 , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);

load resultsBER01.mat;

for i = 1:size(ber,1)
    
    figure
    
    for j = 1:size(ber,2)
        semilogy(SNR,squeeze(ber(i,j,1,:)))
        hold on;
    end
    H = legend('$\mathrm{MI} = 0.05$','$\mathrm{MI} = 0.075$','$\mathrm{MI} = 0.1$');
    set(H,'location','SouthWest');
    set(H,'interpreter','latex')
    
    xlabel('SNR [dB]','interpreter','latex');
    ylim([1e-8 1e0])
    xlim([0 30]);
    
    set(gca,'xtick',SNR);
        
    
    formatFig( gcf ,['.' filesep 'figs' filesep '2017-07-12' filesep 'berFF'  num2str(i)],'en' , figProp );
    
end


fileVector = 2:4;



for k = 1:length(fileVector)
    load(['resultsBER0' num2str(fileVector(k)) '.mat']);
    for i = 1:size(ber,1)

        figure

        for j = 1:size(ber,2)
            semilogy(SNR,squeeze(ber(i,j,1,:)))
            hold on;
        end
        H = legend('$\mathrm{MI} = 0.05$','$\mathrm{MI} = 0.075$','$\mathrm{MI} = 0.1$');
        set(H,'location','SouthWest');
        set(H,'interpreter','latex')

        xlabel('SNR [dB]','interpreter','latex');
        ylim([1e-7 1e0])
        xlim([0 30]);

        set(gca,'xtick',SNR);


    %     formatFig( gcf ,['.' filesep 'figs' filesep '2017-07-12' filesep 'berFF'  num2str(i)],'en' , figProp );

    end
    
end






