clear;
clc;
close all;


addpath(['.' filesep 'results']);

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;

figProp = struct( 'size' , 24 , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);


fileVector = 1:4;


load resultsBER03.mat;

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
%     ylim([1e-6 1e0])
    xlim([0 30]);
    
    set(gca,'xtick',SNR);
        
    
%     formatFig( gcf ,['.' filesep 'figs' filesep '2017-07-12' filesep 'berFF'  num2str(i)],'en' , figProp );
    
end




