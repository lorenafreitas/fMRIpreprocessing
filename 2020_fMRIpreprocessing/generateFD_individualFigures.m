FDdir = '/INSERT HERE THE PARTH TO WHERE YOU SAVED THE FILE allFDs.mat /FramewiseDisplacement/';
allFDsFile = 'allFDs.mat';
M_all = nan(1,30); I_all = nan(1,30); V_all = nan(1,30);

% Select only the tasks you are going to work on:
tasks           = {'RealityFiltering1', 'RealityFiltering2', 'Rest'};
threshold = 0.5; % in mm, as per Power et al. 2014



allFDdata =  dir([FDdir  '*.mat']);
allFDValues = nan(310,length(allFDdata));

for i = 1:length(allFDdata)
    
    load(char(strcat(FDdir,allFDdata(i).name)));
    thresholdLine = ones(length(FD),1)*threshold;
    
    nHighMovFrames = sum(FD>threshold);
    percentHighMov = nHighMovFrames*100/length(FD);
    
    [f2r, percentF2R] = frames2remove((FD>threshold)');
    
    % Plot
     a = figure; plot(FD); hold on;
     plot(thresholdLine, 'LineStyle', ':', 'Linewidth', 3);
     text(1,(threshold + 0.15), 'Threshold', 'FontSize', 15, 'Color', 'r');
     msg = [num2str(percentHighMov,3) '% FD > Thresh.'];
     msg2 = [num2str(percentF2R,3) '% Frames removed'];
     text(50,2.2, msg,'fontSize', 45 , 'Color', [0.4 0.4 0.8]);
     text(50,1.5, msg2,'fontSize', 40 , 'Color', [0.8 0.4 0.4]);
     str = strsplit(allFDdata(i).name,'_');
     title(['FD ' str{2} ' task, Group ' str{5}(1) ', Subj ' str{3}], 'FontSize', 15);
     set(gca, 'Fontsize', 15);
     xlabel('Frames', 'FontSize', 20); ylabel('Framewise Displacement', 'FontSize', 20);
     set(gcf,'color','white'); box off;
     ylim([0 3]);
     
     
     export_fig (['/INSERT HERE THE PATH TO WHERE YOU WANT TO SAVE THE FIGURES/' allFDdata(i).name ], '-pdf');

     close(a);
    
end

close all